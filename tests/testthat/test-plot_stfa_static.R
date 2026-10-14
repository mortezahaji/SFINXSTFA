library(testthat)

# Tests for plot_stfa_static()
# Uses a pre-built minimal STFA_fit result to avoid re-fitting in every test.

# ---- shared fixture --------------------------------------------------------
make_mini_fit <- function(conditions = c("A", "B", "C"), grid_n = 10L) {
  lon  <- seq(0, 100, length.out = grid_n)
  lat  <- seq(0, 100, length.out = grid_n)

  beta_surfaces      <- list()
  significance_masks <- list()
  lower_ci           <- list()
  upper_ci           <- list()

  set.seed(1L)
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

# ---- return type ------------------------------------------------------------

test_that("plot_stfa_static() returns a ggplot object", {
  fit <- make_mini_fit()
  p   <- plot_stfa_static(fit)
  expect_true(inherits(p, "ggplot"))
})

# ---- facet count -----------------------------------------------------------

test_that("plot_stfa_static() produces correct number of facets", {
  fit <- make_mini_fit(conditions = c("Baseline", "Week4", "Response"))
  p   <- plot_stfa_static(fit)

  # Extract facet data from the plot's layout
  built <- ggplot2::ggplot_build(p)
  n_panels <- length(unique(built$layout$layout$PANEL))
  expect_equal(n_panels, 3L)
})

# ---- title -----------------------------------------------------------------

test_that("plot_stfa_static() uses the supplied title", {
  fit   <- make_mini_fit()
  title <- "My Custom Title"
  p     <- plot_stfa_static(fit, title = title)

  expect_equal(p$labels$title, title)
})

# ---- ggsave compatibility --------------------------------------------------

test_that("plot_stfa_static() saves to PNG without error", {
  fit      <- make_mini_fit()
  p        <- plot_stfa_static(fit)
  tmp_file <- tempfile(fileext = ".png")

  expect_no_error(
    ggplot2::ggsave(tmp_file, plot = p, width = 12, height = 4, dpi = 72)
  )
  expect_true(file.exists(tmp_file))
  expect_gt(file.size(tmp_file), 0L)

  unlink(tmp_file)
})

# ---- significance overlay --------------------------------------------------

test_that("plot_stfa_static() sig_overlay = FALSE removes white points layer", {
  fit    <- make_mini_fit()
  p_with <- plot_stfa_static(fit, sig_overlay = TRUE)
  p_no   <- plot_stfa_static(fit, sig_overlay = FALSE)

  # With overlay should have more layers than without
  expect_gt(length(p_with$layers), length(p_no$layers))
})

# ---- NULL images -----------------------------------------------------------

test_that("plot_stfa_static() accepts NULL images without error", {
  fit <- make_mini_fit()
  expect_no_error(plot_stfa_static(fit, images = NULL))
})

# ---- ncol parameter --------------------------------------------------------

test_that("plot_stfa_static() respects ncol parameter", {
  fit <- make_mini_fit(conditions = c("A", "B", "C", "D"))
  p   <- plot_stfa_static(fit, ncol = 2L)

  # Should not error and should still be a ggplot
  expect_true(inherits(p, "ggplot"))
})

# ---- axis labels -----------------------------------------------------------

test_that("plot_stfa_static() has Longitude/Latitude axis labels", {
  fit <- make_mini_fit()
  p   <- plot_stfa_static(fit)

  expect_equal(p$labels$x, "Longitude")
  expect_equal(p$labels$y, "Latitude")
})

# ---- error on bad input ----------------------------------------------------

test_that("plot_stfa_static() errors with missing stfa_result elements", {
  bad_result <- list(conditions = c("A", "B"))  # incomplete
  expect_error(
    plot_stfa_static(bad_result),
    regexp = "missing elements"
  )
})
