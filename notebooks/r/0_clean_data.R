# =====================================================================
# 0_clean_data.R  -  make ONE cleaned dataset that every script uses
# Uses Hannah's updated dataset (10.9.2026), which already has her
# group corrections and no duplicate samples. Stages 1-3 only.
# =====================================================================

library(openxlsx)

RAW <- "data/raw/Cytokines_1236_10.9.2026.xlsx"    # Hannah's file - put it here
OUT <- "data/cleaned/stages1-3_clean.csv"
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

d <- read.xlsx(RAW)

# match the column names the other scripts use
names(d)[names(d) == "participant_id"] <- "ParticipantID"
names(d)[names(d) == "classification"] <- "Dx"
names(d)[names(d) == "stage"]          <- "Stage"
names(d)[names(d) == "batch"]          <- "Batch"
names(d)[names(d) == "IL12P70"]        <- "IL12p70"

CYTO <- c("IFNy","IL1a","IL1b","IL12p70","IL17A","TNFa","IL23",
          "IL4","IL5","IL10","IL13","IL2","IL6","IL8","IL15","TNFb")

d <- d[d$Stage %in% 1:3 & d$Dx %in% c("Control","CFS-C","CFS-R"), ]
d <- d[complete.cases(d[, CYTO]), ]

# check: every participant is in exactly one group
g_per_id <- tapply(d$Dx, d$ParticipantID, function(x) length(unique(x)))
if (any(g_per_id > 1)) {
  cat("\nWARNING - participants listed in more than one group:\n")
  print(names(g_per_id)[g_per_id > 1])
} else cat("\nGroup check: every participant is in exactly one group (no overlap).\n")

# check: one row per participant per stage
key  <- paste(d$ParticipantID, d$Stage)
dups <- unique(key[duplicated(key)])
if (length(dups) > 0) {
  cat("\nWARNING - samples listed more than once:", dups, "\n")
  d <- d[!duplicated(key), ]
} else cat("Duplicate check: one record per participant per stage.\n")

# n table
lab <- c(Control = "Recovered", "CFS-C" = "Moderate", "CFS-R" = "Severe")
cat("\nCleaned n by group and stage:\n")
print(table(factor(lab[d$Dx], levels = lab),
            factor(d$Stage, levels = 1:3, labels = c("Baseline","Mono","6-Month"))))
all3 <- names(which(tapply(d$Stage, d$ParticipantID, function(s) all(1:3 %in% s))))
cat("\nParticipants with all 3 stages:\n")
print(table(factor(lab[d$Dx[match(all3, d$ParticipantID)]], levels = lab)))

keep <- c("ParticipantID", "Stage", "Dx", "Batch", CYTO)
write.csv(d[, keep], OUT, row.names = FALSE)
cat("\nSaved cleaned data:", normalizePath(OUT), "\n")