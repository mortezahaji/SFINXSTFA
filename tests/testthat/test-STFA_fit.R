library(testthat)

# Tests for STFA_fit()
# All tests use make_st_synthetic() — no Seurat required.

# ---- shared fixture --------------------------------------------------------
make_test_data <- function(n = 100L, conditions = c("A", "B", "C"), seed = 42L) {
  make_st_synthetic(
    n_spots    = n,
    conditions = conditions,
    pathways   = c("Response_PW", "Predictor_PW"),
    seed       = seed
  )
}

# ---- basic output structure ------------------------------------------------

test_that("STFA_fit() returns list with required elements", {
  df  <- make_test_data()
  fit <- STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
                  k = 6L, grid_n = 10L)

  expected <- c("model", "conditions", "longitude", "latitude",
                "beta_surfaces", "lower_ci", "upper_ci", "significance_masks")
  expect_true(all(expected %in% names(fit)))
})

test_that("STFA_fit() conditions match input factor levels", {
  df  <- make_test_data(conditions = c("Baseline", "Week4", "Response"))
  fit <- STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
                  k = 6L, grid_n = 10L)

  expect_equal(fit$conditions, c("Baseline", "Week4", "Response"))
})

# ---- grid dimensions -------------------------------------------------------

test_that("STFA_fit() beta surfaces have correct dimensions", {
  df     <- make_test_data()
  grid_n <- 15L
  fit    <- STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
                     k = 6L, grid_n = grid_n)

  for (cond in fit$conditions) {
    expect_equal(dim(fit$beta_surfaces[[cond]]), c(grid_n, grid_n),
                 info = paste("Condition:", cond))
    expect_equal(dim(fit$lower_ci[[cond]]),      c(grid_n, grid_n))
    expect_equal(dim(fit$upper_ci[[cond]]),       c(grid_n, grid_n))
    expect_equal(dim(fit$significance_masks[[cond]]), c(grid_n, grid_n))
  }
})

test_that("STFA_fit() produces one surface per condition", {
  df  <- make_test_data(conditions = c("T1", "T2", "T3", "T4"))
  fit <- STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
                  k = 6L, grid_n = 8L)

  expect_equal(length(fit$beta_surfaces), 4L)
  expect_equal(names(fit$beta_surfaces), c("T1", "T2", "T3", "T4"))
})

# ---- significance mask validity --------------------------------------------

test_that("STFA_fit() significance masks are binary (0/1)", {
  df  <- make_test_data()
  fit <- STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
                  k = 6L, grid_n = 10L)

  for (cond in fit$conditions) {
    vals <- as.integer(fit$significance_masks[[cond]])
    expect_true(all(vals %in% c(0L, 1L)),
                info = paste("Non-binary values in mask for condition:", cond))
  }
})

# ---- CI consistency --------------------------------------------------------

test_that("STFA_fit() lower CI <= beta <= upper CI at each grid cell", {
  df  <- make_test_data(n = 120L)
  fit <- STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
                  k = 6L, grid_n = 10L)

  for (cond in fit$conditions) {
    beta  <- as.numeric(fit$beta_surfaces[[cond]])
    lower <- as.numeric(fit$lower_ci[[cond]])
    upper <- as.numeric(fit$upper_ci[[cond]])

    expect_true(all(lower <= beta  + 1e-9), info = paste(cond, ": lower > beta"))
    expect_true(all(beta  <= upper + 1e-9), info = paste(cond, ": beta > upper"))
  }
})

# ---- grid vectors ----------------------------------------------------------

test_that("STFA_fit() longitude/latitude vectors have length grid_n", {
  df     <- make_test_data()
  grid_n <- 12L
  fit    <- STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
                     k = 6L, grid_n = grid_n)

  expect_equal(length(fit$longitude), grid_n)
  expect_equal(length(fit$latitude),  grid_n)
})

test_that("STFA_fit() grid covers the coordinate range of the input data", {
  df  <- make_test_data()
  fit <- STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
                  k = 6L, grid_n = 10L)

  expect_lte(fit$longitude[1],              min(df$Longitude) + 1e-9)
  expect_gte(fit$longitude[length(fit$longitude)], max(df$Longitude) - 1e-9)
  expect_lte(fit$latitude[1],               min(df$Latitude)  + 1e-9)
  expect_gte(fit$latitude[length(fit$latitude)],  max(df$Latitude)  - 1e-9)
})

# ---- model object ----------------------------------------------------------

test_that("STFA_fit() model element is a bam/gam object", {
  df  <- make_test_data()
  fit <- STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
                  k = 6L, grid_n = 8L)

  expect_true(inherits(fit$model, "gam"),
              info = "model should be a gam/bam object")
})

# ---- error paths -----------------------------------------------------------

test_that("STFA_fit() errors on non-factor Condition", {
  df <- make_test_data()
  df$Condition <- as.character(df$Condition)  # break the factor

  expect_error(
    STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
             k = 6L, grid_n = 8L),
    regexp = "factor"
  )
})

test_that("STFA_fit() errors on missing required columns", {
  df <- make_test_data()
  df$Longitude <- NULL  # remove required column

  expect_error(
    STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
             k = 6L, grid_n = 8L),
    regexp = "missing required columns"
  )
})

test_that("STFA_fit() errors with fewer than 2 conditions", {
  df <- make_test_data(conditions = c("OnlyOne"))

  expect_error(
    STFA_fit(df, Response = "Response_PW", Predictor = "Predictor_PW",
             k = 6L, grid_n = 8L),
    regexp = "at least 2 condition"
  )
})
