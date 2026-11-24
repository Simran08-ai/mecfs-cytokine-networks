# =====================================================================
# make_figures.R  -  one clear figure for each test we ran
#
#   fig1_strength_across_visits.png
#        strength at baseline / mono / 6-month, one line per group,
#        95% CI bars. Top row = all participants, bottom row = only
#        participants with all 3 visits (sensitivity). Global/Pro/Anti/Dual.
#   fig2_curvature_by_group.png

#        curvature D = (baseline + 6-month)/2 - mono for each group, 95% CI
#        (D < 0 = peak at mono, D > 0 = dip at mono). Both analyses.
#   fig3_curvature_differences.png
#        the curvature test itself: Moderate - Recovered and
#        Severe - Recovered, all participants vs all 3 visits, with
#        joint test p and Holm p.
#   fig4_group_differences.png
#        the 36 bootstrap differences (baseline difference, change
#        baseline->mono, change mono->6-month), BH p in the label.
#
# Same data, same strength definition and the same bootstrap as
# bootstrap_differences.R and curvature_test.R (seed 42, 2000 resamples,
# participants resampled within group and keep all their visits), so the
# numbers in the figures are the same numbers as in those CSVs.
# Saves in figures/final_figures/
# =====================================================================

set.seed(42)
NBOOT <- 2000
d <- read.csv("data/cleaned/stages1-3_clean.csv", check.names = FALSE, stringsAsFactors = FALSE)
OUT <- "figures/final_figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CYTO <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23",
          "IL4","IL5","IL10","IL13","IL2","IL6","IL8","IL15","TNFb")
PRO  <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23")
ANTI <- c("IL4","IL5","IL10","IL13")
DUAL <- c("IL2","IL6","IL8","IL15","TNFb")
MEAS <- c("Global","Pro","Anti","Dual")
MEAS_LAB <- c(Global = "Global (all 16)", Pro = "Pro-inflammatory", Anti = "Anti-inflammatory", Dual = "Dual")
GROUPS <- c(Recovered = "Control", Moderate = "CFS-C", Severe = "CFS-R")
GCOL   <- c(Recovered = "#2E7D32", Moderate = "#E67E22", Severe = "#8E44AD")
GPCH   <- c(Recovered = 19, Moderate = 15, Severe = 17)
STAGE_LAB <- c("Baseline", "Mono", "6-Month")

summaries <- function(m) {
  cm <- suppressWarnings(cor(m, method = "spearman")); diag(cm) <- 0
  s  <- rowSums(abs(cm))                       # all |r|, no threshold
  c(Global = mean(s), Pro = mean(s[PRO]), Anti = mean(s[ANTI]), Dual = mean(s[DUAL]))
}
q  <- function(x, p) quantile(x, p, na.rm = TRUE, names = FALSE)
pb <- function(x) { x <- x[!is.na(x)]; min(1, 2 * min(mean(x <= 0), mean(x >= 0))) }   # two-sided bootstrap p
curv <- function(x) (x[, 1] + x[, 3]) / 2 - x[, 2]

# ---- bootstrap: same loop order as curvature_test.R ----
run <- function(dd) {
  lookup <- list(); ids <- list()
  for (g in names(GROUPS)) {
    dg <- dd[dd$Dx == GROUPS[[g]], ]
    ids[[g]] <- unique(dg$ParticipantID)
    lookup[[g]] <- lapply(1:3, function(s) {
      m <- as.matrix(dg[dg$Stage == s, CYTO]); rownames(m) <- dg$ParticipantID[dg$Stage == s]; m })
  }
  gs <- function(g, draw) sapply(1:3, function(s) {
    m <- lookup[[g]][[s]]; summaries(m[draw[draw %in% rownames(m)], , drop = FALSE]) })
  obs  <- lapply(names(GROUPS), function(g) gs(g, ids[[g]])); names(obs) <- names(GROUPS)
  boot <- lapply(names(GROUPS), function(g) array(NA, c(NBOOT, 4, 3), dimnames = list(NULL, MEAS, NULL)))
  names(boot) <- names(GROUPS)
  for (b in 1:NBOOT) for (g in names(GROUPS)) {
    draw <- sample(ids[[g]], length(ids[[g]]), replace = TRUE)
    boot[[g]][b, , ] <- gs(g, draw)
  }
  n <- sapply(names(GROUPS), function(g) sapply(lookup[[g]], nrow))   # 3 stages x 3 groups
  list(obs = obs, boot = boot, n = n)
}
cat("Bootstrap, all participants...\n");  A <- run(d)
all3 <- names(which(tapply(d$Stage, d$ParticipantID, function(s) all(1:3 %in% s))))
cat("Bootstrap, all 3 visits only...\n"); B <- run(d[d$ParticipantID %in% all3, ])
AN <- list("All participants" = A, "All 3 visits only" = B)

# ---- tables behind the figures ----
traj <- do.call(rbind, lapply(names(AN), function(a) do.call(rbind, lapply(names(GROUPS), function(g)
  do.call(rbind, lapply(MEAS, function(m) data.frame(Analysis = a, Group = g, Measure = m, Stage = STAGE_LAB,
    n = AN[[a]]$n[, g], Strength = AN[[a]]$obs[[g]][m, ],
    CI_low  = apply(AN[[a]]$boot[[g]][, m, ], 2, q, .025),
    CI_high = apply(AN[[a]]$boot[[g]][, m, ], 2, q, .975))))))))

cv <- do.call(rbind, lapply(names(AN), function(a) do.call(rbind, lapply(names(GROUPS), function(g) {
  bc <- t(apply(AN[[a]]$boot[[g]], 1, curv))          # NBOOT x 4
  data.frame(Analysis = a, Group = g, Measure = MEAS, D = curv(AN[[a]]$obs[[g]]),
             CI_low = apply(bc, 2, q, .025), CI_high = apply(bc, 2, q, .975)) }))))

cd <- list(); jt <- list()
for (a in names(AN)) for (m in MEAS) {
  X <- AN[[a]]
  bc <- sapply(names(GROUPS), function(g) apply(X$boot[[g]][, m, ], 1, function(r) (r[1] + r[3]) / 2 - r[2]))
  oc <- sapply(names(GROUPS), function(g) curv(X$obs[[g]])[m])
  dObs <- c(oc[2] - oc[1], oc[3] - oc[1])
  dB <- cbind(bc[, 2] - bc[, 1], bc[, 3] - bc[, 1]); dB <- dB[complete.cases(dB), , drop = FALSE]
  S <- cov(dB); Si <- solve(S); W <- as.numeric(t(dObs) %*% Si %*% dObs)
  cen <- sweep(dB, 2, colMeans(dB)); Wb <- rowSums((cen %*% Si) * cen)
  jt[[length(jt) + 1]] <- data.frame(Analysis = a, Measure = m, Wald = W, p_joint = mean(Wb >= W))
  p <- c(pb(dB[, 1]), pb(dB[, 2]))
  cd[[length(cd) + 1]] <- data.frame(Analysis = a, Measure = m,
    Comparison = c("Moderate - Recovered", "Severe - Recovered"), Diff = dObs,
    CI_low = apply(dB, 2, q, .025), CI_high = apply(dB, 2, q, .975), p_boot = p, p_Holm = p.adjust(p, "holm"))
}
cd <- do.call(rbind, cd); jt <- do.call(rbind, jt)

TYPES <- c(base = "Baseline difference", ch12 = "Change baseline -> mono", ch23 = "Change mono -> 6-month")
PAIRS <- list(c("Moderate", "Recovered"), c("Severe", "Recovered"), c("Severe", "Moderate"))
stat <- function(x, type, m) {        # x: 4x3 obs or NBOOT x 4 x 3 boot
  if (length(dim(x)) == 2) switch(type, base = x[m, 1], ch12 = x[m, 2] - x[m, 1], ch23 = x[m, 3] - x[m, 2])
  else switch(type, base = x[, m, 1], ch12 = x[, m, 2] - x[, m, 1], ch23 = x[, m, 3] - x[, m, 2])
}
gd <- do.call(rbind, lapply(names(TYPES), function(t) do.call(rbind, lapply(PAIRS, function(pr)
  do.call(rbind, lapply(MEAS, function(m) {
    bd <- stat(A$boot[[pr[1]]], t, m) - stat(A$boot[[pr[2]]], t, m)
    data.frame(Contrast = TYPES[[t]], Comparison = paste(pr[1], "-", pr[2]), Measure = m,
               Estimate = stat(A$obs[[pr[1]]], t, m) - stat(A$obs[[pr[2]]], t, m),
               CI_low = q(bd, .025), CI_high = q(bd, .975), p_boot = pb(bd)) }))))))
gd$p_BH <- p.adjust(gd$p_boot, "BH")

r2 <- function(df) { i <- sapply(df, is.numeric); df[i] <- lapply(df[i], round, 3); df }
write.csv(r2(traj), file.path(OUT, "fig1_numbers.csv"), row.names = FALSE)
write.csv(r2(cv),   file.path(OUT, "fig2_numbers.csv"), row.names = FALSE)
write.csv(merge(r2(cd), r2(jt)), file.path(OUT, "fig3_numbers.csv"), row.names = FALSE)
write.csv(r2(gd),   file.path(OUT, "fig4_numbers.csv"), row.names = FALSE)

fp <- function(p) ifelse(p < .001, "<.001", sub("^0", "", sprintf("%.3f", round(p, 3))))
ngroup <- function(X, g) paste(X$n[, g], collapse = "/")

# =====================================================================
# Figure 1 - strength across visits
# =====================================================================
png(file.path(OUT, "fig1_strength_across_visits.png"), width = 4200, height = 2300, res = 220)
par(mfrow = c(2, 4), mar = c(3, 4.2, 3, 0.8), oma = c(3.5, 2.5, 3, 0))
for (a in names(AN)) for (m in MEAS) {
  tt <- traj[traj$Measure == m, ]
  yl <- range(c(tt$CI_low, tt$CI_high)) + c(-0.2, 0.2)
  tt <- tt[tt$Analysis == a, ]
  plot(NA, xlim = c(0.7, 3.3), ylim = yl, xaxt = "n", xlab = "", ylab = "strength (sum |r|)", main = MEAS_LAB[m])
  axis(1, at = 1:3, labels = STAGE_LAB); grid(nx = NA, ny = NULL, col = "gray88")
  off <- c(Recovered = -0.08, Moderate = 0, Severe = 0.08)
  for (g in names(GROUPS)) {
    s <- tt[tt$Group == g, ]; x <- 1:3 + off[g]
    arrows(x, s$CI_low, x, s$CI_high, angle = 90, code = 3, length = 0.04, col = GCOL[g], lwd = 1.3)
    lines(x, s$Strength, col = GCOL[g], lwd = 2.5); points(x, s$Strength, pch = GPCH[g], col = GCOL[g], cex = 1.4)
  }
  if (m == "Global") mtext(a, side = 2, line = 4.8, font = 2, cex = 0.95)
}
mtext("Network strength across the three visits (point = estimate, bar = 95% bootstrap CI)", outer = TRUE, font = 2, cex = 1.2, line = 1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
legend("bottom", ncol = 3, bty = "n", lwd = 2.5, pch = GPCH, col = GCOL, cex = 1.05, text.width = 0.28,
       legend = paste0(names(GROUPS), " (n = ", sapply(names(GROUPS), function(g) ngroup(A, g)),
                       " all; ", sapply(names(GROUPS), function(g) ngroup(B, g)), " all 3 visits)"))
dev.off()

# =====================================================================
# Figure 2 - curvature for each group
# =====================================================================
png(file.path(OUT, "fig2_curvature_by_group.png"), width = 4200, height = 2000, res = 220)
par(mfrow = c(2, 4), mar = c(4, 7.5, 3, 1), oma = c(2.5, 2.5, 3, 0))
xl <- range(c(cv$CI_low, cv$CI_high)) + c(-0.2, 0.2)
for (a in names(AN)) for (m in MEAS) {
  s <- cv[cv$Analysis == a & cv$Measure == m, ]; y <- 3:1
  plot(NA, xlim = xl, ylim = c(0.5, 3.5), yaxt = "n", ylab = "", xlab = "D", main = MEAS_LAB[m])
  axis(2, at = y, labels = s$Group, las = 1); abline(v = 0, lty = 2, col = "gray50")
  grid(nx = NULL, ny = NA, col = "gray90")
  segments(s$CI_low, y, s$CI_high, y, col = GCOL[s$Group], lwd = 2.5)
  points(s$D, y, pch = GPCH[s$Group], col = GCOL[s$Group], cex = 1.6)
  if (m == "Global") mtext(a, side = 2, line = 6.3, font = 2, cex = 0.95)
}
mtext("Curvature D = (baseline + 6-month)/2 - mono   (D < 0: peak at mono,  D > 0: dip at mono)",
      outer = TRUE, font = 2, cex = 1.2, line = 1)
mtext("point = estimate, line = 95% bootstrap CI; if the line crosses the dashed 0 line, no clear dip or peak",
      side = 1, outer = TRUE, line = 1, cex = 0.85)
dev.off()

# =====================================================================
# Figure 3 - curvature test (group differences in curvature)
# =====================================================================
ACOL <- c("All participants" = "black", "All 3 visits only" = "#2E86C1")
png(file.path(OUT, "fig3_curvature_differences.png"), width = 4200, height = 1500, res = 220)
par(mfrow = c(1, 4), mar = c(4.5, 10.5, 4.5, 1), oma = c(2.5, 0, 2.5, 0))
xl <- range(c(cd$CI_low, cd$CI_high)) + c(-0.3, 0.3)
for (m in MEAS) {
  s <- cd[cd$Measure == m, ]
  s <- s[order(s$Comparison, s$Analysis != "All participants"), ]
  y <- c(4.2, 3.5, 1.7, 1.0)
  plot(NA, xlim = xl, ylim = c(0.4, 4.8), yaxt = "n", ylab = "", xlab = "difference in D",
       main = paste0(MEAS_LAB[m], "\njoint p: ", fp(jt$p_joint[jt$Measure == m & jt$Analysis == "All participants"]),
                     " (all) / ", fp(jt$p_joint[jt$Measure == m & jt$Analysis == "All 3 visits only"]), " (all 3)"),
       cex.main = 1)
  axis(2, at = c(3.85, 1.35), labels = c("Moderate -\nRecovered", "Severe -\nRecovered"), las = 1, tick = FALSE)
  abline(v = 0, lty = 2, col = "gray50"); grid(nx = NULL, ny = NA, col = "gray90")
  segments(s$CI_low, y, s$CI_high, y, col = ACOL[s$Analysis], lwd = 2.5)
  points(s$Diff, y, pch = 19, col = ACOL[s$Analysis], cex = 1.5)
  text(xl[2], y + 0.28, paste0("Holm p = ", fp(s$p_Holm)), adj = 1, cex = 0.8, col = ACOL[s$Analysis])
}
mtext("Curvature test: does the mono dip/peak differ from Recovered?  (point = difference, line = 95% bootstrap CI)",
      outer = TRUE, font = 2, cex = 1.15, line = 0.8)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
legend("bottom", horiz = TRUE, bty = "n", lwd = 2.5, pch = 19, col = ACOL, legend = names(ACOL), cex = 1.05)
dev.off()

# =====================================================================
# Figure 4 - the 36 bootstrap group differences (all participants)
# =====================================================================
PCOL <- c("Moderate - Recovered" = "#E67E22", "Severe - Recovered" = "#8E44AD", "Severe - Moderate" = "#555555")
png(file.path(OUT, "fig4_group_differences.png"), width = 4200, height = 2300, res = 220)
par(mfrow = c(1, 3), mar = c(4.5, 11, 3.5, 1), oma = c(3, 0, 3, 0))
xl <- range(c(gd$CI_low, gd$CI_high)) + c(-0.3, 2.2)
for (t in TYPES) {
  s <- gd[gd$Contrast == t, ]; y <- rev(seq_len(nrow(s))) + rep(c(0.6, 0.3, 0), each = 4)
  plot(NA, xlim = xl, ylim = c(0.5, max(y) + 0.5), yaxt = "n", ylab = "", xlab = "difference", main = t)
  axis(2, at = y, labels = paste0(s$Measure, "  "), las = 1, tick = FALSE, cex.axis = 0.9)
  for (k in 1:3) mtext(names(PCOL)[k], side = 2, at = mean(y[(4 * k - 3):(4 * k)]), line = 4.6, cex = 0.7, font = 2, col = PCOL[k])
  abline(v = 0, lty = 2, col = "gray50"); grid(nx = NULL, ny = NA, col = "gray90")
  sig <- s$CI_low > 0 | s$CI_high < 0
  segments(s$CI_low, y, s$CI_high, y, col = PCOL[s$Comparison], lwd = 2.2)
  points(s$Estimate, y, pch = ifelse(sig, 19, 1), col = PCOL[s$Comparison], cex = 1.5, lwd = 2)
  text(xl[2], y, paste0("BH p ", fp(s$p_BH)), adj = 1, cex = 0.72, col = "gray30")
}
mtext("Bootstrap group differences in strength (all participants; filled point = 95% CI does not cross 0)",
      outer = TRUE, font = 2, cex = 1.15, line = 1)
mtext("Change = later visit minus earlier visit, compared between groups. BH = corrected across all 36 tests.",
      side = 1, outer = TRUE, line = 1.2, cex = 0.85)
dev.off()

cat("\nJoint curvature test:\n"); print(r2(jt), row.names = FALSE)
cat("\nCurvature differences:\n"); print(r2(cd), row.names = FALSE)
cat("\nDone. 4 figures + 4 number tables saved in:", normalizePath(OUT), "\n")
