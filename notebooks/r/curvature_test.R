# =====================================================================
# curvature_test.R  -  the "mono dip vs peak" test Lenny asked for
#
# For each group:  Curvature = (baseline + 6-month) / 2  -  mono
#   positive = network strength DIPS at mono
#   negative = network strength PEAKS at mono
#
# Strength definition (stated explicitly):
#   For each group at each stage, Spearman correlations between all 16
#   cytokines. Each cytokine's strength = sum of the ABSOLUTE correlations
#   with the other 15 cytokines, using ALL correlations (no |r| > .75
#   cut-off - the .75 cut-off is only used to draw the network graphs).
#   Global = mean strength of all 16; Pro / Anti / Dual = mean of the
#   cytokines in that class.
#
# Bootstrap:
#   2000 resamples. Within each group, participants are drawn with
#   replacement; each drawn participant keeps all of their observations
#   from stages 1-3. Networks are rebuilt at each stage from the resample.
#   The curvature and the group differences are computed inside every
#   replicate, so the three visits stay linked.
#
# Tests (primary = Global strength):
#   1) Joint test of the two curvature differences
#      (Moderate - Recovered, Severe - Recovered) together:
#      Wald statistic W = d' S^-1 d, S = bootstrap covariance of d.
#      p from the bootstrap (share of centred replicates with W_b >= W)
#      and from chi-square with 2 df for comparison.
#   2) The two individual contrasts: 95% percentile CI and two-sided
#      bootstrap p = 2 * min( P(diff <= 0), P(diff >= 0) ), Holm-adjusted.
#   Pro / Anti / Dual are reported the same way as exploratory.
#
# Sensitivity: everything repeated using only participants who have all
# three stages.
# =====================================================================

set.seed(42)
NBOOT <- 2000
DATA  <- "data/cleaned/stages1-3_clean.csv"     # made by 0_clean_data.R
OUT   <- "reports/curvature"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(DATA, check.names = FALSE, stringsAsFactors = FALSE)

CYTO <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23",
          "IL4","IL5","IL10","IL13","IL2","IL6","IL8","IL15","TNFb")
PRO  <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23")
ANTI <- c("IL4","IL5","IL10","IL13")
DUAL <- c("IL2","IL6","IL8","IL15","TNFb")
MEAS <- c("Global","Pro","Anti","Dual")
GROUPS <- c(Recovered = "Control", Moderate = "CFS-C", Severe = "CFS-R")

summaries <- function(m) {
  cm <- suppressWarnings(cor(m, method = "spearman")); diag(cm) <- 0
  s  <- rowSums(abs(cm))                     # all |r|, no threshold
  c(Global = mean(s), Pro = mean(s[PRO]), Anti = mean(s[ANTI]), Dual = mean(s[DUAL]))
}

run_analysis <- function(dd, label) {
  cat("\n=====================================================================\n")
  cat(" ", label, "\n")
  cat("=====================================================================\n")

  lookup <- list(); ids <- list()
  for (g in names(GROUPS)) {
    dg <- dd[dd$Dx == GROUPS[[g]], ]
    ids[[g]] <- unique(dg$ParticipantID)
    lookup[[g]] <- lapply(1:3, function(s) {
      m <- as.matrix(dg[dg$Stage == s, CYTO]); rownames(m) <- dg$ParticipantID[dg$Stage == s]; m })
  }
  group_summ <- function(g, draw) sapply(1:3, function(s) {
    m <- lookup[[g]][[s]]; summaries(m[draw[draw %in% rownames(m)], , drop = FALSE]) })
  curv <- function(x) (x[, 1] + x[, 3]) / 2 - x[, 2]          # x: 4 x 3

  cat("\nn (baseline / mono / 6-month):\n")
  for (g in names(GROUPS))
    cat(sprintf("  %-10s %3d / %3d / %3d   (participants: %d)\n", g,
                nrow(lookup[[g]][[1]]), nrow(lookup[[g]][[2]]), nrow(lookup[[g]][[3]]), length(ids[[g]])))

  # observed
  obs_s <- lapply(names(GROUPS), function(g) group_summ(g, ids[[g]])); names(obs_s) <- names(GROUPS)
  obs_c <- sapply(obs_s, curv)                                    # 4 measures x 3 groups

  cat("\nGlobal strength by stage:\n")
  gt <- t(sapply(obs_s, function(x) round(x["Global", ], 2))); colnames(gt) <- c("Baseline","Mono","6-Month")
  print(gt)
  cat("\nCurvature = (baseline + 6-month)/2 - mono   (+ = dip at mono, - = peak at mono):\n")
  print(round(obs_c, 2))

  # bootstrap
  bc <- array(NA, c(NBOOT, 4, 3), dimnames = list(NULL, MEAS, names(GROUPS)))
  for (b in 1:NBOOT) for (g in names(GROUPS)) {
    draw <- sample(ids[[g]], length(ids[[g]]), replace = TRUE)
    bc[b, , g] <- curv(group_summ(g, draw))
  }

  rows <- list(); joint <- list()
  for (m in MEAS) {
    dObs <- c(obs_c[m, "Moderate"] - obs_c[m, "Recovered"], obs_c[m, "Severe"] - obs_c[m, "Recovered"])
    dB   <- cbind(bc[, m, "Moderate"] - bc[, m, "Recovered"], bc[, m, "Severe"] - bc[, m, "Recovered"])
    dB   <- dB[complete.cases(dB), , drop = FALSE]

    # joint Wald test
    S  <- cov(dB); Si <- solve(S)
    W  <- as.numeric(t(dObs) %*% Si %*% dObs)
    cen <- sweep(dB, 2, colMeans(dB))
    Wb <- rowSums((cen %*% Si) * cen)
    joint[[m]] <- data.frame(Analysis = label, Measure = m,
                             Wald = round(W, 2), p_boot = round(mean(Wb >= W), 3),
                             p_chisq_df2 = round(pchisq(W, df = 2, lower.tail = FALSE), 3))

    # individual contrasts
    for (k in 1:2) {
      p <- min(1, 2 * min(mean(dB[, k] <= 0), mean(dB[, k] >= 0)))
      rows[[length(rows) + 1]] <- data.frame(
        Analysis = label, Measure = m,
        Comparison = c("Moderate - Recovered", "Severe - Recovered")[k],
        Curv_diff = round(dObs[k], 2),
        CI_low  = round(quantile(dB[, k], 0.025, names = FALSE), 2),
        CI_high = round(quantile(dB[, k], 0.975, names = FALSE), 2),
        p_boot  = round(p, 3))
    }
  }
  ind <- do.call(rbind, rows)
  ind$p_Holm <- NA
  for (m in MEAS) { i <- ind$Measure == m; ind$p_Holm[i] <- round(p.adjust(ind$p_boot[i], "holm"), 3) }
  jt <- do.call(rbind, joint)

  cat("\nJoint test of both curvature differences (Global = primary):\n")
  print(jt[, -1], row.names = FALSE)
  cat("\nIndividual curvature differences (Holm within each measure):\n")
  print(ind[, -1], row.names = FALSE)

  list(individual = ind, joint = jt,
       curvature = data.frame(Analysis = label, Measure = rownames(obs_c), round(obs_c, 3)))
}

# ---- main analysis: everyone with any of stages 1-3 ----
main <- run_analysis(d, "All participants (stages 1-3)")

# ---- sensitivity: only participants with all three stages ----
all3 <- names(which(tapply(d$Stage, d$ParticipantID, function(s) all(1:3 %in% s))))
sens <- run_analysis(d[d$ParticipantID %in% all3, ], "Sensitivity: participants with all 3 stages")

write.csv(rbind(main$individual, sens$individual), file.path(OUT, "curvature_individual.csv"), row.names = FALSE)
write.csv(rbind(main$joint, sens$joint),           file.path(OUT, "curvature_joint.csv"),      row.names = FALSE)
write.csv(rbind(main$curvature, sens$curvature),   file.path(OUT, "curvature_observed.csv"),   row.names = FALSE)
cat("\nSaved results in:", normalizePath(OUT), "\n")
