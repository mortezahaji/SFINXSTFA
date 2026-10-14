library(testthat)

# Tests for run_spatiotemporal()
# All tests bypass Seurat by using the 'data' argument with synthetic data.

# ---- shared fixture --------------------------------------------------------
make_pipeline_data <- function(n = 150L,
                                conditions = c("Baseline", "Week4", "Response"),
                                seed = 42L) {
  make_st_synthetic(
    n_spots    = n,
    conditions = conditions,
    pathways   = c("HP_LUPUS_NEPHRITIS", "GOBP_NEUTROPHIL"),
    seed       = seed
  )
}

# ---- return type and element names -----------------------------------------

test_that("run_spatiotemporal() returns list with 4 named elements", {
  df     <- make_pipeline_data()
  result <- run_spatiotemporal(
    data      = df,
    response  = "HP_LUPUS_NEPHRITIS",
    predictor = "GOBP_NEUTROPHIL",
    k         = 6L,
    grid_n    = 10L,
    title     = "Test_Pipeline"
  )

  expect_true(is.list(result))
  expect_equal(
    sort(names(result)),
    sort(c("data", "model_fit", "static_plot", "interactive_plot"))
  )
})

# ---- data element ----------------------------------------------------------

test_that("run_spatiotemporal() data element is a data.frame", {
  df     <- make_pipeline_data()
  result <- run_spatiotemporal(
    data      = df,
    response  = "HP_LUPUS_NEPHRITIS",
    predictor = "GOBP_NEUTROPHIL",
    k         = 6L,
    grid_n    = 10L
  )

  expect_true(is.data.frame(result$data))
  expect_equal(nrow(result$data), nrow(df))
})

# ---- model_fit element -----------------------------------------------------

test_that("run_spatiotemporal() model_fit has correct structure", {
  df     <- make_pipeline_data()
  result <- run_spatiotemporal(
    data      = df,
    response  = "HP_LUPUS_NEPHRITIS",
    predictor = "GOBP_NEUTROPHIL",
    k         = 6L,
    grid_n    = 10L
  )

  expected <- c("model", "conditions", "longitude", "latitude",
                "beta_surfaces", "lower_ci", "upper_ci", "significance_masks")
  expect_true(all(expected %in% names(result$model_fit)))
})

# ---- static_plot element ---------------------------------------------------

test_that("run_spatiotemporal() static_plot is a ggplot", {
  df     <- make_pipeline_data()
  result <- run_spatiotemporal(
    data      = df,
    response  = "HP_LUPUS_NEPHRITIS",
    predictor = "GOBP_NEUTROPHIL",
    k         = 6L,
    grid_n    = 10L
  )

  expect_true(inherits(result$static_plot, "ggplot"))
})

# ---- interactive_plot element ----------------------------------------------

test_that("run_spatiotemporal() interactive_plot is a plotly/htmlwidget", {
  df     <- make_pipeline_data()
  result <- run_spatiotemporal(
    data      = df,
    response  = "HP_LUPUS_NEPHRITIS",
    predictor = "GOBP_NEUTROPHIL",
    k         = 6L,
    grid_n    = 10L
  )

  expect_true(
    inherits(result$interactive_plot, "plotly") ||
    inherits(result$interactive_plot, "htmlwidget")
  )
})

# ---- output_dir saves files ------------------------------------------------

test_that("run_spatiotemporal() saves files when output_dir is provided", {
  df      <- make_pipeline_data()
  tmp_dir <- tempfile()

  result <- run_spatiotemporal(
    data       = df,
    response   = "HP_LUPUS_NEPHRITIS",
    predictor  = "GOBP_NEUTROPHIL",
    k          = 6L,
    grid_n     = 10L,
    output_dir = tmp_dir,
    title      = "Test_Output"
  )

  expect_true(dir.exists(tmp_dir))

  # Static PNG
  static_file <- file.path(tmp_dir, "Test_Output_static.png")
  expect_true(file.exists(static_file))
  expect_gt(file.size(static_file), 0L)

  # Interactive 2D HTML
  html_2d <- file.path(tmp_dir, "Test_Output_interactive_2d.html")
  expect_true(file.exists(html_2d))
  expect_gt(file.size(html_2d), 0L)

  # Interactive 3D HTML
  html_3d <- file.path(tmp_dir, "Test_Output_interactive_3d.html")
  expect_true(file.exists(html_3d))
  expect_gt(file.size(html_3d), 0L)

  unlink(tmp_dir, recursive = TRUE)
})

# ---- error: neither seurat_list nor data provided --------------------------

test_that("run_spatiotemporal() errors when neither seurat_list nor data provided", {
  expect_error(
    run_spatiotemporal(response = "R", predictor = "P"),
    regexp = "seurat_list.*data|data.*seurat_list"
  )
})

# ---- error: both seurat_list and data provided -----------------------------

test_that("run_spatiotemporal() errors when both seurat_list and data provided", {
  df <- make_pipeline_data()
  expect_error(
    run_spatiotemporal(seurat_list = list(), data = df,
                       response = "HP_LUPUS_NEPHRITIS", predictor = "GOBP_NEUTROPHIL"),
    regexp = "only one"
  )
})

# ---- error: fewer than 2 conditions ----------------------------------------

test_that("run_spatiotemporal() errors with fewer than 2 conditions", {
  df <- make_pipeline_data(conditions = c("OnlyOne"))

  expect_error(
    run_spatiotemporal(
      data      = df,
      response  = "HP_LUPUS_NEPHRITIS",
      predictor = "GOBP_NEUTROPHIL",
      k         = 6L,
      grid_n    = 8L
    ),
    regexp = "at least 2 condition"
  )
})

# ---- error: missing Condition column ---------------------------------------

test_that("run_spatiotemporal() errors when Condition column is missing", {
  df <- make_pipeline_data()
  df$Condition <- NULL

  expect_error(
    run_spatiotemporal(
      data      = df,
      response  = "HP_LUPUS_NEPHRITIS",
      predictor = "GOBP_NEUTROPHIL",
      k         = 6L,
      grid_n    = 8L
    ),
    regexp = "missing required columns"
  )
})

# ---- error: missing response column ----------------------------------------

test_that("run_spatiotemporal() errors when response column is missing", {
  df <- make_pipeline_data()

  expect_error(
    run_spatiotemporal(
      data      = df,
      response  = "NON_EXISTENT_PATHWAY",
      predictor = "GOBP_NEUTROPHIL",
      k         = 6L,
      grid_n    = 8L
    ),
    regexp = "missing required columns"
  )
})

# ---- conditions in model match input ---------------------------------------

test_that("run_spatiotemporal() model conditions match data conditions", {
  conditions <- c("Pre", "Mid", "Post")
  df <- make_pipeline_data(conditions = conditions)

  result <- run_spatiotemporal(
    data      = df,
    response  = "HP_LUPUS_NEPHRITIS",
    predictor = "GOBP_NEUTROPHIL",
    k         = 6L,
    grid_n    = 10L
  )

  expect_equal(result$model_fit$conditions, conditions)
})
