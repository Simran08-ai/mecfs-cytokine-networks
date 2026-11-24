# =====================================================================
# ME/CFS cytokine networks + weighted edge strength + 95% CIs
# Stages 1-3, three groups (Recovered / Moderate / Severe)
# Saves all 18 figures as PNGs + the results table as a CSV
# =====================================================================

library(igraph)
set.seed(42)   # same confidence intervals every run

# ---- where everything gets saved ----
OUT <- "figures/final_report"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- load raw stages 1-3 ----
# cleaned data from 0_clean_data.R (duplicates handled, stages 1-3, 3 groups)
d <- read.csv("data/cleaned/stages1-3_clean.csv",
              check.names = FALSE, stringsAsFactors = FALSE)

# all 16 cytokines, and which class each one is in
CYTO <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23",
          "IL4","IL5","IL10","IL13","IL2","IL6","IL8","IL15","TNFb")
PRO  <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23")
ANTI <- c("IL4","IL5","IL10","IL13")
DUAL <- c("IL2","IL6","IL8","IL15","TNFb")

CLASS <- c(setNames(rep("pro",length(PRO)), PRO),
           setNames(rep("anti",length(ANTI)), ANTI),
           setNames(rep("dual",length(DUAL)), DUAL))
COLOR <- c(pro="lightcoral", anti="lightblue", dual="lightgray")
PTCOL <- c(pro="#C0392B", anti="#2E86C1", dual="#7F8C8D")

# groups and stages
GROUPS <- list(Recovered = "Control", Moderate = "CFS-C", Severe = "CFS-R")
STAGES <- c(Baseline = 1, Mono = 2, "6-Month" = 3)

# grab one group's cytokine matrix at one stage
get_sub <- function(dx, stage){
  sub <- d[d$Dx == dx & d$Stage == stage, CYTO]
  sub[complete.cases(sub), ]
}

# =====================================================================
# PART 1 - network graphs (the 9 graphs)
# =====================================================================
# edge kept if |r| > 0.75, edge width scaled by correlation strength,
# Fruchterman-Reingold layout (same as Kurt's).
draw_network <- function(dx, stage, title){
  sub <- get_sub(dx, stage)
  cm  <- cor(sub, method = "spearman"); diag(cm) <- 0
  A   <- abs(cm); A[A < 0.75] <- 0              # keep only strong links
  g   <- graph_from_adjacency_matrix(A, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)
  E(g)$width    <- abs(E(g)$weight)^2 * 5        # thicker = stronger
  V(g)$color    <- COLOR[CLASS[V(g)$name]]
  layout_coords <- layout_with_fr(g, niter = 15000)
  plot(g, layout = layout_coords, vertex.size = 20,
       vertex.label.cex = 0.8, vertex.label.color = "black",
       margin = 0,
       main = paste0(title, " (spearman, |r| >0.75) n=", nrow(sub)))
}

# make all 9 and save each one
for (st in names(STAGES))
  for (gn in names(GROUPS)){
    fn <- file.path(OUT, paste0("net_", gn, "_", st, ".png"))
    png(fn, width = 1800, height = 1400, res = 200)
    draw_network(GROUPS[[gn]], STAGES[[st]], paste(gn, st))
    dev.off()
  }

# =====================================================================
# PART 2 - weighted edge strength by class, with 95% CIs
# =====================================================================
# strength = sum of each cytokine's absolute correlations (full network),
# then averaged across the cytokines in a class. CI = bootstrap over people.

class_strength <- function(sub, class){
  cm <- cor(sub, method = "spearman"); diag(cm) <- 0
  s  <- rowSums(abs(cm))                          # per-cytokine strength
  mean(s[names(s) %in% names(CLASS)[CLASS == class]])
}

boot_ci <- function(dx, stage, class, nboot = 2000){
  sub <- get_sub(dx, stage)
  point <- class_strength(sub, class)
  vals <- replicate(nboot, {
    idx <- sample(nrow(sub), nrow(sub), replace = TRUE)
    class_strength(sub[idx, ], class)
  })
  c(strength = point,
    lo = quantile(vals, 0.025, names = FALSE),
    hi = quantile(vals, 0.975, names = FALSE),
    n  = nrow(sub))
}

# build the full table (pro / anti / dual x 3 stages x 3 groups)
rows <- list()
for (class in c("pro","anti","dual"))
  for (st in names(STAGES))
    for (gn in names(GROUPS)){
      r <- boot_ci(GROUPS[[gn]], STAGES[[st]], class)
      rows[[length(rows)+1]] <- data.frame(
        Class = class, Stage = st, Group = gn, n = r["n"],
        Strength = round(r["strength"], 2),
        CI_low = round(r["lo"], 2), CI_high = round(r["hi"], 2))
    }
results <- do.call(rbind, rows)
print(results, row.names = FALSE)
write.csv(results, file.path(OUT, "strength_ci_table.csv"), row.names = FALSE)

# plot: one per class per stage, point + 95% CI error bars, and save each
CLASS_NAME <- c(pro = "Pro-inflammatory", anti = "Anti-inflammatory", dual = "Dual")
for (class in c("pro","anti","dual"))
  for (st in names(STAGES)){
    sub <- results[results$Class == class & results$Stage == st, ]
    fn  <- file.path(OUT, paste0("ci_", class, "_", st, ".png"))
    png(fn, width = 1500, height = 1000, res = 200)
    par(mar = c(5, 5, 4, 2))
    plot(1:3, sub$Strength,
         ylim = c(min(sub$CI_low) - 0.5, max(sub$CI_high) + 0.5),
         xlim = c(0.5, 3.5), xaxt = "n", pch = 19, cex = 1.6,
         col = PTCOL[class], xlab = "", ylab = "strength (sum |r|)",
         main = paste0(CLASS_NAME[class], " strength - ", st, "\n(point + 95% CI)"))
    arrows(1:3, sub$CI_low, 1:3, sub$CI_high, angle = 90, code = 3,
           length = 0.08, lwd = 2, col = PTCOL[class])
    axis(1, at = 1:3, labels = paste0(sub$Group, "\n(n=", sub$n, ")"),
         padj = 0.5)
    grid(nx = NA, ny = NULL, col = "gray85")
    dev.off()
  }

cat("\nDone. 18 PNGs + strength_ci_table.csv saved in:", normalizePath(OUT), "\n")
