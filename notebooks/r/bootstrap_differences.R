# =====================================================================
# Bootstrap differences between groups - for Lenny
#
# Three groups kept separate: Recovered (Control), Moderate (CFS-C),
# Severe (CFS-R). Stages 1-3 only (baseline, mono, 6-month).
# Nobody needs stage 6 data to be included.
#
# Principal network summaries (weighted degree strength, |r| network):
#   Global - mean strength over all 16 cytokines
#   Pro    - mean strength of the 7 pro-inflammatory cytokines
#   Anti   - mean strength of the 4 anti-inflammatory cytokines
#   Dual   - mean strength of the 5 dual cytokines
#
# Contrasts:
#   A) baseline differences between groups
#   B) change baseline -> mono, compared between groups
#   C) change mono -> 6-month, compared between groups
#
# Resampling: participants are resampled WITHIN each group, and each
# selected participant keeps all their observations across stages 1-3
# (linked / participant-level bootstrap).
#
# Output: printed table + CSV with estimate, 95% CI, bootstrap p-value,
# and a Benjamini-Hochberg adjusted p-value.
# =====================================================================

set.seed(42)
NBOOT <- 2000

OUT <- "reports/bootstrap_differences"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- data ----
# cleaned data from 0_clean_data.R (duplicates handled, stages 1-3, 3 groups)
d <- read.csv("data/cleaned/stages1-3_clean.csv",
              check.names = FALSE, stringsAsFactors = FALSE)

CYTO <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23",
          "IL4","IL5","IL10","IL13","IL2","IL6","IL8","IL15","TNFb")
PRO  <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23")
ANTI <- c("IL4","IL5","IL10","IL13")
DUAL <- c("IL2","IL6","IL8","IL15","TNFb")

GROUPS <- c(Recovered = "Control", Moderate = "CFS-C", Severe = "CFS-R")
STAGES <- 1:3
STAGE_LAB <- c("Baseline", "Mono", "6-Month")

d <- d[d$Stage %in% STAGES & d$Dx %in% GROUPS, ]
d <- d[complete.cases(d[, CYTO]), ]
d <- d[!duplicated(d[, c("ParticipantID", "Stage")]), ]   # one row per person per stage

# lookup: for each group and stage, a matrix of cytokines with rownames = ParticipantID
lookup <- list()
ids    <- list()
for (g in names(GROUPS)) {
  dg <- d[d$Dx == GROUPS[[g]], ]
  ids[[g]] <- unique(dg$ParticipantID)          # everyone in the group with any of stages 1-3
  lookup[[g]] <- lapply(STAGES, function(s) {
    m <- as.matrix(dg[dg$Stage == s, CYTO])
    rownames(m) <- dg$ParticipantID[dg$Stage == s]
    m
  })
}

# ---- the four principal summaries from one cytokine matrix ----
summaries <- function(m) {
  cm <- suppressWarnings(cor(m, method = "spearman")); diag(cm) <- 0
  s  <- rowSums(abs(cm))
  c(Global = mean(s), Pro = mean(s[PRO]), Anti = mean(s[ANTI]), Dual = mean(s[DUAL]))
}

# summaries for one group at all 3 stages, given a set of participant IDs
# (an ID drawn twice contributes its row twice at every stage it has)
group_summaries <- function(g, draw) {
  sapply(STAGES, function(s) {
    m    <- lookup[[g]][[s]]
    keep <- draw[draw %in% rownames(m)]
    summaries(m[keep, , drop = FALSE])
  })   # 4 x 3 matrix: measures x stages
}

# ---- observed values ----
obs <- lapply(names(GROUPS), function(g) group_summaries(g, ids[[g]]))
names(obs) <- names(GROUPS)

# ---- bootstrap ----
boot <- lapply(names(GROUPS), function(g) array(NA, c(NBOOT, 4, 3)))
names(boot) <- names(GROUPS)
for (b in 1:NBOOT) {
  for (g in names(GROUPS)) {
    draw <- sample(ids[[g]], length(ids[[g]]), replace = TRUE)
    boot[[g]][b, , ] <- group_summaries(g, draw)
  }
}

MEAS <- c("Global", "Pro", "Anti", "Dual")

# value of a quantity for a group: level at a stage, or change between two stages
stat <- function(x, type, m) {      # x: 4x3 (obs) or n x 4 x 3 (boot)
  if (length(dim(x)) == 2) {
    switch(type, base = x[m, 1], ch12 = x[m, 2] - x[m, 1], ch23 = x[m, 3] - x[m, 2])
  } else {
    switch(type, base = x[, m, 1], ch12 = x[, m, 2] - x[, m, 1], ch23 = x[, m, 3] - x[, m, 2])
  }
}
TYPE_LAB <- c(base = "Baseline difference",
              ch12 = "Change Baseline->Mono, difference",
              ch23 = "Change Mono->6-Month, difference")

PAIRS <- list(c("Moderate", "Recovered"), c("Severe", "Recovered"), c("Severe", "Moderate"))

rows <- list()
for (type in names(TYPE_LAB)) for (pr in PAIRS) for (mi in seq_along(MEAS)) {
  a <- pr[1]; r <- pr[2]
  est <- stat(obs[[a]], type, mi) - stat(obs[[r]], type, mi)
  bd  <- stat(boot[[a]], type, mi) - stat(boot[[r]], type, mi)
  bd  <- bd[!is.na(bd)]
  p   <- min(1, 2 * min(mean(bd <= 0), mean(bd >= 0)))   # two-sided bootstrap p
  rows[[length(rows) + 1]] <- data.frame(
    Contrast = TYPE_LAB[[type]],
    Comparison = paste(a, "-", r),
    Measure = MEAS[mi],
    Estimate = round(est, 2),
    CI_low = round(quantile(bd, 0.025, names = FALSE), 2),
    CI_high = round(quantile(bd, 0.975, names = FALSE), 2),
    p_boot = round(p, 3))
}
res <- do.call(rbind, rows)
res$p_BH <- round(p.adjust(res$p_boot, method = "BH"), 3)
res$CI_excludes_0 <- ifelse(res$CI_low > 0 | res$CI_high < 0, "yes", "no")

# ---- sample sizes ----
cat("\nParticipants per group (any of stages 1-3) and per stage:\n")
for (g in names(GROUPS))
  cat(sprintf("  %-10s total=%3d   baseline=%3d  mono=%3d  6-month=%3d\n", g,
              length(ids[[g]]), nrow(lookup[[g]][[1]]), nrow(lookup[[g]][[2]]), nrow(lookup[[g]][[3]])))

cat("\nObserved strength (Global / Pro / Anti / Dual) by stage:\n")
for (g in names(GROUPS)) {
  o <- round(obs[[g]], 2); colnames(o) <- STAGE_LAB
  cat("\n", g, "\n"); print(o)
}

cat("\nBootstrap differences (", NBOOT, " resamples, linked by participant):\n\n", sep = "")
print(res, row.names = FALSE)

write.csv(res, file.path(OUT, "bootstrap_differences.csv"), row.names = FALSE)
cat("\nSaved:", normalizePath(file.path(OUT, "bootstrap_differences.csv")), "\n")
