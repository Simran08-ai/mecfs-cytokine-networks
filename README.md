# ME/CFS Cytokine Networks

R code for the ME/CFS cytokine correlation network analysis at DePaul University, with Dr. Leonard A. Jason and Dr. Jacob Furst. This work extends the earlier cytokine network analysis and uses Hannah's corrected dataset.

For each group (Recovered, Moderate, Severe) at each visit (baseline, mono, 6 months), I build a Spearman correlation network of 16 cytokines and look at how the network changes across the three visits.

**The data is not included in this repo.** It is participant data from the study and is kept private.

## What the scripts do

Run them in this order, from the main project folder (the folder that has `data/` and `notebooks/` in it).

| # | Script | What it does |
|---|--------|--------------|
| 1 | `0_clean_data.R` | Reads Hannah's file and makes one cleaned dataset (stages 1-3, 3 groups, one record per person per visit) |
| 2 | `mecfs_networks_and_strength.R` | Draws the 9 network graphs (\|r\| > 0.75) and the Pro / Anti / Dual strength graphs with 95% CIs |
| 3 | `bootstrap_differences.R` | Bootstrap group differences (36 comparisons, BH correction) |
| 4 | `curvature_test.R` | Curvature test for the dip/peak at mono (joint Wald test, Holm correction), also with only people who have all 3 visits |
| 5 | `make_figures.R` | Makes the 4 main figures: strength across visits, curvature by group, curvature test, group differences |
| 6 | `density_clustering.R` | Network density, clustering and weighted clustering, with bootstrap CIs |
| 7 | `clustering_check.R` | Checks how cytokines with fewer than 2 links are counted, and how often the bootstrap gives zero clustering |

## Methods in short

- **Strength:** each cytokine's strength is the sum of its absolute correlations with the other 15 (all correlations, no cut-off). Global = average of all 16; Pro (7), Anti (4) and Dual (5) = average within each class.
- **The 0.75 cut-off** is only used for the network graphs, density and unweighted clustering. Strength and weighted clustering use all correlations.
- **Bootstrap:** 2000 resamples, participants drawn with replacement within each group, seed 42.
- **Curvature:** D = (baseline + 6-month) / 2 - mono. Negative = peak at mono, positive = dip at mono.

## How to run

1. Install R and the packages: `install.packages(c("igraph", "openxlsx"))`
2. Put the data file at `data/raw/Cytokines_1236_10.9.2026.xlsx` (not included here)
3. Set the working directory to the main project folder and run the scripts in the order above, for example:
   ```r
   source("notebooks/r/0_clean_data.R")
   ```

Results and figures are saved in `figures/` and `reports/`.

## Author

Simran Chauhan, DePaul University
