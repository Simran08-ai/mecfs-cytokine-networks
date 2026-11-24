# =====================================================================
# density_clustering.R  -  network density and clustering
# Same setup as the weighted edge strength analysis:
#   - same cleaned data (made by 0_clean_data.R)
#   - 3 groups x 3 stages, Spearman correlations of the 16 cytokines
#   - Global (all 16) and Pro / Anti / Dual (average of the cytokines in that class)
#   - bootstrap within each group and stage, 2000 resamples, 95% percentile CI
#
# Measures
#   1. Density (|r| > 0.75 network):
#        each cytokine: number of strong links / 15
#        Global = average over 16 cytokines = strong links / 120
#   2. Clustering (|r| > 0.75 network):
#        each cytokine: of its linked partners, the share that are also linked
#        to each other (igraph local transitivity; fewer than 2 links = 0)
#   3. Weighted clustering (all correlations, no threshold):
#        each cytokine: average of (|r_ij| * |r_ik| * |r_jk|)^(1/3) over all
#        pairs of its partners j, k  (Onnela weighted clustering; |r| is
#        already between 0 and 1)
#
# Saves results in figures/density_clustering/
# =====================================================================

library(igraph)
set.seed(42)
NBOOT <- 2000
THR   <- 0.75

d <- read.csv("data/cleaned/stages1-3_clean.csv", check.names = FALSE, stringsAsFactors = FALSE)
OUT <- "figures/density_clustering"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CYTO <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23",
          "IL4","IL5","IL10","IL13","IL2","IL6","IL8","IL15","TNFb")
PRO  <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23")
ANTI <- c("IL4","IL5","IL10","IL13")
DUAL <- c("IL2","IL6","IL8","IL15","TNFb")
GROUPS <- c(Recovered = "Control", Moderate = "CFS-C", Severe = "CFS-R")
STAGES <- c(Baseline = 1, Mono = 2, "6-Month" = 3)
GCOL   <- c(Recovered = "#2E7D32", Moderate = "#E67E22", Severe = "#8E44AD")

# per-cytokine values -> Global / Pro / Anti / Dual
by_class <- function(v) c(Global = mean(v), Pro = mean(v[PRO]), Anti = mean(v[ANTI]), Dual = mean(v[DUAL]))

# all three measures from one set of samples
metrics <- function(m) {
  cm <- suppressWarnings(cor(m, method = "spearman")); diag(cm) <- 0
  cm[is.na(cm)] <- 0
  W <- abs(cm)
  A <- (W > THR) * 1                                   # strong-link network (0/1)
  dens <- rowSums(A) / 15
  g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  clus <- transitivity(g, type = "local", isolates = "zero"); names(clus) <- CYTO
  clus[is.na(clus)] <- 0
  Wc <- W^(1/3)
  wclus <- diag(Wc %*% Wc %*% Wc) / (15 * 14)          # Onnela, complete weighted network
  c(Density = by_class(dens), Clustering = by_class(clus), WeightedClustering = by_class(wclus))
}

rows <- list()
for (st in names(STAGES)) for (g in names(GROUPS)) {
  sub <- as.matrix(d[d$Dx == GROUPS[[g]] & d$Stage == STAGES[[st]], CYTO])
  est <- metrics(sub)
  bt  <- replicate(NBOOT, metrics(sub[sample(nrow(sub), nrow(sub), replace = TRUE), ]))
  lo  <- apply(bt, 1, quantile, 0.025, na.rm = TRUE)
  hi  <- apply(bt, 1, quantile, 0.975, na.rm = TRUE)
  for (k in names(est)) {
    p <- strsplit(k, ".", fixed = TRUE)[[1]]
    rows[[length(rows) + 1]] <- data.frame(Measure = p[1], Class = p[2], Stage = st, Group = g,
      n = nrow(sub), Value = round(est[k], 3), CI_low = round(lo[k], 3), CI_high = round(hi[k], 3))
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(OUT, "density_clustering_table.csv"), row.names = FALSE)

cat("\nGlobal values (all 16 cytokines):\n")
print(res[res$Class == "Global", c("Measure", "Stage", "Group", "n", "Value", "CI_low", "CI_high")], row.names = FALSE)

# one figure per measure: rows = Global/Pro/Anti/Dual, columns = stages
TITLE <- c(Density = "Network density (|r| > 0.75)", Clustering = "Clustering (|r| > 0.75)",
           WeightedClustering = "Weighted clustering (all correlations)")
for (ms in names(TITLE)) {
  png(file.path(OUT, paste0(ms, ".png")), width = 3000, height = 3200, res = 250)
  par(mfrow = c(4, 3), mar = c(3, 4.2, 2.5, 0.8), oma = c(1, 0, 3, 0))
  r <- res[res$Measure == ms, ]
  for (cl in c("Global", "Pro", "Anti", "Dual")) {
    yl <- range(c(r$CI_low[r$Class == cl], r$CI_high[r$Class == cl]), na.rm = TRUE)
    for (st in names(STAGES)) {
      s <- r[r$Class == cl & r$Stage == st, ]
      plot(1:3, s$Value, ylim = yl, xlim = c(0.5, 3.5), xaxt = "n", pch = 19, cex = 1.5,
           col = GCOL[s$Group], xlab = "", ylab = cl, main = paste(cl, "-", st), cex.main = 1)
      arrows(1:3, s$CI_low, 1:3, s$CI_high, angle = 90, code = 3, length = 0.05, lwd = 2, col = GCOL[s$Group])
      axis(1, at = 1:3, labels = paste0(s$Group, "\n(n=", s$n, ")"), padj = 0.4, cex.axis = 0.8)
      grid(nx = NA, ny = NULL, col = "gray88")
    }
  }
  mtext(paste0(TITLE[ms], "  -  point = value, bar = 95% bootstrap CI"), outer = TRUE, font = 2, cex = 1.1, line = 1)
  dev.off()
}
cat("\nDone. Table and 3 figures saved in:", normalizePath(OUT), "\n")
