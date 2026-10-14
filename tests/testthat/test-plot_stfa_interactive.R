library(testthat)

# Tests for plot_stfa_interactive()
# Uses the same minimal fixture as plot_stfa_static tests.

# ---- shared fixture --------------------------------------------------------
make_mini_fit <- function(conditions = c("A", "B", "C"), grid_n = 10L) {
  lon <- seq(0, 100, length.out = grid_n)
  lat <- seq(0, 100, length.out = grid_n)

  beta_surfaces      <- list()
  significance_masks <- list()
  lower_ci           <- list()
  upper_ci           <- list()

  set.seed(7L)
  for (cond in conditions) {
    bm <- matrix(stats::rnorm(grid_n * grid_n), nrow = grid_n)
    se <- matrix(abs(stats::rnorm(grid_n * grid_n, sd = 0.3)), nrow = grid_n)
    beta_surfaces[[cond]]      <- bm
    lower_ci[[cond]]           <- bm - 1.96 * se
    upper_ci[[cond]]           <- bm + 1.96 * se
    significance_masks[[cond]] <- matrix(
      as.integer((bm - 1.96 * se > 0) | (bm + 1.96 * se < 0)),
      nrow = grid_n
    )
  }

  list(
    conditions         = conditions,
    longitude          = lon,
    latitude           = lat,
    beta_surfaces      = beta_surfaces,
    lower_ci           = lower_ci,
    upper_ci           = upper_ci,
    significance_masks = significance_masks
  )
}

# ---- return type: 2D -------------------------------------------------------

test_that("plot_stfa_interactive() mode='2d' returns a plotly/htmlwidget", {
  fit <- make_mini_fit()
  p   <- plot_stfa_interactive(fit, mode = "2d")

  expect_true(inherits(p, "plotly") || inherits(p, "htmlwidget"),
              info = "Should return a plotly or htmlwidget object")
})

# ---- return type: 3D -------------------------------------------------------

test_that("plot_stfa_interactive() mode='3d' returns a plotly/htmlwidget", {
  fit <- make_mini_fit()
  p   <- plot_stfa_interactive(fit, mode = "3d")

  expect_true(inherits(p, "plotly") || inherits(p, "htmlwidget"))
})

# ---- slider steps equal number of conditions: 2D ---------------------------

test_that("plot_stfa_interactive() 2D slider has n_conditions steps", {
  conditions <- c("Baseline", "Week4", "Response")
  fit <- make_mini_fit(conditions = conditions)
  p   <- plot_stfa_interactive(fit, mode = "2d")

  # Extract slider steps from plotly layout
  n_steps <- length(p$x$layout$sliders[[1]]$steps)
  expect_equal(n_steps, length(conditions))
})

# ---- slider steps equal number of conditions: 3D ---------------------------

test_that("plot_stfa_interactive() 3D slider has n_conditions steps", {
  conditions <- c("T0", "T1", "T2", "T3")
  fit <- make_mini_fit(conditions = conditions)
  p   <- plot_stfa_interactive(fit, mode = "3d")

  n_steps <- length(p$x$layout$sliders[[1]]$steps)
  expect_equal(n_steps, length(conditions))
})

# ---- slider step labels match condition names ------------------------------

test_that("plot_stfa_interactive() slider labels match condition names", {
  conditions <- c("Baseline", "Week4", "Response")
  fit <- make_mini_fit(conditions = conditions)
  p   <- plot_stfa_interactive(fit, mode = "2d")

  labels <- vapply(
    p$x$layout$sliders[[1]]$steps,
    function(s) s$label,
    character(1L)
  )
  expect_equal(labels, conditions)
})

# ---- trace count: 2D -------------------------------------------------------

test_that("plot_stfa_interactive() 2D has 2 traces per condition (heatmap + contour)", {
  conditions <- c("A", "B", "C")
  fit <- make_mini_fit(conditions = conditions)
  p   <- plot_stfa_interactive(fit, mode = "2d")

  n_traces <- length(p$x$data)
  expect_equal(n_traces, length(conditions) * 2L)
})

# ---- trace count: 3D -------------------------------------------------------

test_that("plot_stfa_interactive() 3D has 1 surface trace per condition", {
  conditions <- c("A", "B")
  fit <- make_mini_fit(conditions = conditions)
  p   <- plot_stfa_interactive(fit, mode = "3d")

  n_traces <- length(p$x$data)
  expect_equal(n_traces, length(conditions))
})

# ---- title is set ----------------------------------------------------------

test_that("plot_stfa_interactive() respects title parameter", {
  fit   <- make_mini_fit()
  mytitle <- "Custom Explorer Title"
  p     <- plot_stfa_interactive(fit, title = mytitle)

  plot_title <- p$x$layout$title$text
  expect_equal(plot_title, mytitle)
})

# ---- mode default is "2d" --------------------------------------------------

test_that("plot_stfa_interactive() defaults to mode='2d'", {
  fit <- make_mini_fit()
  p   <- plot_stfa_interactive(fit)  # no mode argument

  # 2D mode produces 2 traces per condition (heatmap + contour)
  n_traces <- length(p$x$data)
  expect_equal(n_traces, length(fit$conditions) * 2L)
})

# ---- error on bad input ----------------------------------------------------

test_that("plot_stfa_interactive() errors with missing stfa_result elements", {
  bad_result <- list(conditions = c("A", "B"))
  expect_error(
    plot_stfa_interactive(bad_result),
    regexp = "missing elements"
  )
})

# ---- invalid mode ----------------------------------------------------------

test_that("plot_stfa_interactive() errors on invalid mode", {
  fit <- make_mini_fit()
  expect_error(
    plot_stfa_interactive(fit, mode = "4d"),
    regexp = "should be one of"
  )
})
