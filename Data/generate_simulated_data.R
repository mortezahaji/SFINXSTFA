## =============================================================================
## generate_simulated_data.R
## Generates all simulated datasets for testing SFinxSDFA v0.2.0 functions.
## Run this script once to populate the Data/ folder with CSV files.
##
## Datasets produced:
##   1. st_3condition_450spots.csv   — main dataset (3 conditions x 150 spots)
##   2. st_2condition_200spots.csv   — minimal 2-condition dataset
##   3. st_5condition_750spots.csv   — stress-test: 5 conditions x 150 spots
##   4. st_singlespot_grid.csv       — regular 10x10 grid layout per condition
##   5. st_multipathway.csv          — 5 pathways for binary operator testing
## =============================================================================

set.seed(42)

## ---------------------------------------------------------------------------
## Helper: generate spatially structured pathway scores
## Produces scores with a sinusoidal spatial gradient that shifts across
## conditions, giving STFA_fit() real signal to detect.
## ---------------------------------------------------------------------------
make_scores <- function(lon, lat, condition_idx, pathway_idx, sd = 1.2) {
  n <- length(lon)
  gradient <- sin(lon / 100 * pi) * cos(lat / 100 * pi) * condition_idx
  base <- rnorm(n, mean = 0, sd = sd)
  base + gradient * (0.5 + 0.3 * pathway_idx)
}

## ============================================================================
## Dataset 1: st_3condition_450spots.csv
## Primary test dataset — used to exercise all five new functions.
##   - 3 conditions: Baseline, Week4, Response
##   - 150 spots per condition (450 rows total)
##   - 2 pathways: GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION, HP_LUPUS_NEPHRITIS
##   - Columns: Longitude, Latitude, Condition, PatientID + 2 pathway columns
## ============================================================================
cat("Generating Dataset 1: st_3condition_450spots.csv ...\n")

conditions_3 <- c("Baseline", "Week4", "Response")
n_spots      <- 150L

rows_3 <- lapply(seq_along(conditions_3), function(i) {
  cond <- conditions_3[i]
  lon  <- runif(n_spots, 0, 100)
  lat  <- runif(n_spots, 0, 100)
  data.frame(
    Longitude                              = round(lon, 4),
    Latitude                               = round(lat, 4),
    Condition                              = cond,
    PatientID                              = paste0("Patient_", i),
    GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION = round(make_scores(lon, lat, i, 1), 4),
    HP_LUPUS_NEPHRITIS                     = round(make_scores(lon, lat, i, 2), 4),
    stringsAsFactors                       = FALSE
  )
})

df1 <- do.call(rbind, rows_3)
df1$Condition <- factor(df1$Condition, levels = conditions_3)
rownames(df1) <- NULL

write.csv(df1, file = "Data/st_3condition_450spots.csv", row.names = FALSE)
cat("  Saved: Data/st_3condition_450spots.csv  [", nrow(df1), "rows x",
    ncol(df1), "cols]\n")

## ============================================================================
## Dataset 2: st_2condition_200spots.csv
## Minimal 2-condition dataset — for quick unit tests and edge-case checks.
##   - 2 conditions: Pre, Post
##   - 100 spots per condition (200 rows total)
##   - 2 pathways
## ============================================================================
cat("Generating Dataset 2: st_2condition_200spots.csv ...\n")

conditions_2 <- c("Pre", "Post")
n_spots_2    <- 100L

rows_2 <- lapply(seq_along(conditions_2), function(i) {
  cond <- conditions_2[i]
  lon  <- runif(n_spots_2, 0, 100)
  lat  <- runif(n_spots_2, 0, 100)
  data.frame(
    Longitude                              = round(lon, 4),
    Latitude                               = round(lat, 4),
    Condition                              = cond,
    PatientID                              = paste0("Patient_", i),
    GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION = round(make_scores(lon, lat, i, 1), 4),
    HP_LUPUS_NEPHRITIS                     = round(make_scores(lon, lat, i, 2), 4),
    stringsAsFactors                       = FALSE
  )
})

df2 <- do.call(rbind, rows_2)
df2$Condition <- factor(df2$Condition, levels = conditions_2)
rownames(df2) <- NULL

write.csv(df2, file = "Data/st_2condition_200spots.csv", row.names = FALSE)
cat("  Saved: Data/st_2condition_200spots.csv  [", nrow(df2), "rows x",
    ncol(df2), "cols]\n")

## ============================================================================
## Dataset 3: st_5condition_750spots.csv
## Stress-test dataset — 5 conditions for testing slider with many steps.
##   - 5 conditions: T0, T1, T2, T3, T4 (time series)
##   - 150 spots per condition (750 rows total)
##   - 2 pathways
##   - Multiple patients (2 per condition for random effect testing)
## ============================================================================
cat("Generating Dataset 3: st_5condition_750spots.csv ...\n")

conditions_5 <- c("T0", "T1", "T2", "T3", "T4")
n_spots_5    <- 150L

rows_5 <- lapply(seq_along(conditions_5), function(i) {
  cond     <- conditions_5[i]
  lon      <- runif(n_spots_5, 0, 100)
  lat      <- runif(n_spots_5, 0, 100)
  # Two patients per condition (first 75 = PatientA, last 75 = PatientB)
  patient  <- c(rep(paste0("PatientA_", i), 75), rep(paste0("PatientB_", i), 75))
  data.frame(
    Longitude                              = round(lon, 4),
    Latitude                               = round(lat, 4),
    Condition                              = cond,
    PatientID                              = patient,
    GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION = round(make_scores(lon, lat, i, 1), 4),
    HP_LUPUS_NEPHRITIS                     = round(make_scores(lon, lat, i, 2, sd = 1.0), 4),
    stringsAsFactors                       = FALSE
  )
})

df3 <- do.call(rbind, rows_5)
df3$Condition <- factor(df3$Condition, levels = conditions_5)
rownames(df3) <- NULL

write.csv(df3, file = "Data/st_5condition_750spots.csv", row.names = FALSE)
cat("  Saved: Data/st_5condition_750spots.csv  [", nrow(df3), "rows x",
    ncol(df3), "cols]\n")

## ============================================================================
## Dataset 4: st_grid_layout.csv
## Regular grid layout — tests SFinx_UO and SFinx_BO on evenly spaced spots.
## Mimics a Visium capture array layout (10x10 grid per condition).
##   - 3 conditions: Baseline, Week4, Response
##   - 100 spots per condition on a 10x10 regular grid (300 rows total)
##   - 2 pathways with clear spatial structure
## ============================================================================
cat("Generating Dataset 4: st_grid_layout.csv ...\n")

grid_seq <- seq(5, 95, length.out = 10)
grid_xy  <- expand.grid(Longitude = grid_seq, Latitude = grid_seq)

rows_grid <- lapply(seq_along(conditions_3), function(i) {
  cond <- conditions_3[i]
  lon  <- grid_xy$Longitude
  lat  <- grid_xy$Latitude
  data.frame(
    Longitude                              = round(lon, 4),
    Latitude                               = round(lat, 4),
    Condition                              = cond,
    PatientID                              = paste0("Patient_", i),
    GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION = round(make_scores(lon, lat, i, 1, sd = 0.8), 4),
    HP_LUPUS_NEPHRITIS                     = round(make_scores(lon, lat, i, 2, sd = 0.8), 4),
    stringsAsFactors                       = FALSE
  )
})

df4 <- do.call(rbind, rows_grid)
df4$Condition <- factor(df4$Condition, levels = conditions_3)
rownames(df4) <- NULL

write.csv(df4, file = "Data/st_grid_layout.csv", row.names = FALSE)
cat("  Saved: Data/st_grid_layout.csv  [", nrow(df4), "rows x",
    ncol(df4), "cols]\n")

## ============================================================================
## Dataset 5: st_multipathway.csv
## Multi-pathway dataset — for testing SFinx_BO binary operators (Union,
## Intersection, Difference, Dominance) across pairs of pathways.
##   - 3 conditions: Baseline, Week4, Response
##   - 150 spots per condition (450 rows total)
##   - 5 pathways with different spatial patterns and activation rates
## ============================================================================
cat("Generating Dataset 5: st_multipathway.csv ...\n")

pathways_5 <- c(
  "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION",
  "HP_LUPUS_NEPHRITIS",
  "GOBP_INFLAMMATORY_RESPONSE",
  "HP_PROTEINURIA",
  "GOBP_CYTOKINE_PRODUCTION"
)

# Pathway-specific SD and mean offsets to create varied activation patterns
pw_params <- list(
  list(sd = 1.5, mean_shift = 0.0),   # PW1: broad, centred
  list(sd = 1.2, mean_shift = 0.2),   # PW2: moderate
  list(sd = 1.8, mean_shift = -0.3),  # PW3: high variance, left-shifted
  list(sd = 0.9, mean_shift = 0.5),   # PW4: narrow, right-shifted (high activation)
  list(sd = 1.3, mean_shift = 0.0)    # PW5: moderate
)

rows_mp <- lapply(seq_along(conditions_3), function(i) {
  cond <- conditions_3[i]
  lon  <- runif(n_spots, 0, 100)
  lat  <- runif(n_spots, 0, 100)

  pw_scores <- lapply(seq_along(pathways_5), function(j) {
    p     <- pw_params[[j]]
    n     <- n_spots
    grad  <- sin(lon / 100 * pi) * cos(lat / 100 * pi) * i
    score <- rnorm(n, mean = p$mean_shift, sd = p$sd) + grad * (0.4 + 0.2 * j)
    round(score, 4)
  })
  names(pw_scores) <- pathways_5

  df_row <- data.frame(
    Longitude = round(lon, 4),
    Latitude  = round(lat, 4),
    Condition = cond,
    PatientID = paste0("Patient_", i),
    stringsAsFactors = FALSE
  )
  cbind(df_row, as.data.frame(pw_scores))
})

df5 <- do.call(rbind, rows_mp)
df5$Condition <- factor(df5$Condition, levels = conditions_3)
rownames(df5) <- NULL

write.csv(df5, file = "Data/st_multipathway.csv", row.names = FALSE)
cat("  Saved: Data/st_multipathway.csv  [", nrow(df5), "rows x",
    ncol(df5), "cols]\n")

## ============================================================================
## Summary
## ============================================================================
cat("\n=== All datasets saved to Data/ ===\n")
cat(sprintf("  %-40s  %4d rows  %d cols\n",
    "st_3condition_450spots.csv",  nrow(df1), ncol(df1)))
cat(sprintf("  %-40s  %4d rows  %d cols\n",
    "st_2condition_200spots.csv",  nrow(df2), ncol(df2)))
cat(sprintf("  %-40s  %4d rows  %d cols\n",
    "st_5condition_750spots.csv",  nrow(df3), ncol(df3)))
cat(sprintf("  %-40s  %4d rows  %d cols\n",
    "st_grid_layout.csv",          nrow(df4), ncol(df4)))
cat(sprintf("  %-40s  %4d rows  %d cols\n",
    "st_multipathway.csv",         nrow(df5), ncol(df5)))
cat("\nColumn reference:\n")
cat("  Longitude, Latitude     — spot coordinates in [0, 100]\n")
cat("  Condition               — factor; levels define time points\n")
cat("  PatientID               — donor ID for random effects in STFA_fit()\n")
cat("  [pathway columns]       — z-scored pathway activity (mean ~0, SD ~1.2)\n")
