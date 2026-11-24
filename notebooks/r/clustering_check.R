# =====================================================================
# clustering_check.R  -  two checks Dr. Jason asked for on thresholded clustering
#
# 1) Are cytokines with fewer than 2 links counted as 0 when averaging?
#    Prints, for Severe 6-month, each cytokine's number of strong links and its
#    clustering value, and shows that Global = sum of all 16 values / 16.
# 2) How often do the bootstrap networks give zero clustering?
#    Reruns the same bootstrap as density_clustering.R (same seed 42, same
#    order, 2000 resamples) so the CIs come out the same, and counts how many
#    of the 2000 bootstrap networks have clustering = 0, for every group and visit.
#
# Uses the cleaned file from 0_clean_data.R
# Saves figures/density_clustering/clustering_zero_check.csv
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

by_class <- function(v) c(Global = mean(v), Pro = mean(v[PRO]), Anti = mean(v[ANTI]), Dual = mean(v[DUAL]))

# exactly the same function as in density_clustering.R
metrics <- function(m) {
  cm <- suppressWarnings(cor(m, method = "spearman")); diag(cm) <- 0
  cm[is.na(cm)] <- 0
  W <- abs(cm)
  A <- (W > THR) * 1
  dens <- rowSums(A) / 15
  g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  clus <- transitivity(g, type = "local", isolates = "zero"); names(clus) <- CYTO
  clus[is.na(clus)] <- 0
  Wc <- W^(1/3)
  wclus <- diag(Wc %*% Wc %*% Wc) / (15 * 14)
  c(Density = by_class(dens), Clustering = by_class(clus), WeightedClustering = by_class(wclus))
}

# ---- check 1: how cytokines with fewer than 2 links are handled ----
sev6 <- as.matrix(d[d$Dx == "CFS-R" & d$Stage == 3, CYTO])
cm <- cor(sev6, method = "spearman"); diag(cm) <- 0
A  <- (abs(cm) > THR) * 1
g  <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
cl <- transitivity(g, type = "local", isolates = "zero"); cl[is.na(cl)] <- 0
cat("\nCheck 1 - Severe 6-month, observed network\n")
print(data.frame(Cytokine = CYTO, StrongLinks = rowSums(A), Clustering = round(cl, 3)), row.names = FALSE)
cat("Strong links in the network:", sum(A) / 2, "\n")
cat("Cytokines with fewer than 2 links:", sum(rowSums(A) < 2), "(all counted as 0)\n")
cat("Global clustering = sum of all 16 values / 16 =", round(sum(cl) / 16, 3), "\n")

# ---- check 2: share of bootstrap networks with zero clustering ----
rows <- list()
for (st in names(STAGES)) for (gname in names(GROUPS)) {
  sub <- as.matrix(d[d$Dx == GROUPS[[gname]] & d$Stage == STAGES[[st]], CYTO])
  est <- metrics(sub)
  bt  <- replicate(NBOOT, metrics(sub[sample(nrow(sub), nrow(sub), replace = TRUE), ]))
  for (cl_name in c("Global", "Pro", "Anti", "Dual")) {
    k <- paste0("Clustering.", cl_name)
    rows[[length(rows) + 1]] <- data.frame(Stage = st, Group = gname, Class = cl_name, n = nrow(sub),
      Observed = round(est[k], 3),
      CI_low = round(quantile(bt[k, ], 0.025, names = FALSE), 3),
      CI_high = round(quantile(bt[k, ], 0.975, names = FALSE), 3),
      Boot_zero = sum(bt[k, ] == 0),
      Boot_zero_pct = round(100 * mean(bt[k, ] == 0), 1))
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(OUT, "clustering_zero_check.csv"), row.names = FALSE)

cat("\nCheck 2 - Global clustering: how many of the 2000 bootstrap networks had clustering = 0\n")
print(res[res$Class == "Global", ], row.names = FALSE)
cat("\nSevere 6-month, all classes:\n")
print(res[res$Group == "Severe" & res$Stage == "6-Month", ], row.names = FALSE)
cat("\nSaved:", normalizePath(file.path(OUT, "clustering_zero_check.csv")), "\n")
