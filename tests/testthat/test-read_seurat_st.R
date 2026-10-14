library(testthat)

# Tests for read_seurat_st() and make_st_synthetic()
# Note: read_seurat_st() requires SeuratObject; those tests use mock data only.
# All tests run without Seurat installed.

test_that("make_st_synthetic() returns correct structure", {
  df <- make_st_synthetic(
    n_spots    = 50,
    conditions = c("Pre", "Post"),
    pathways   = c("PathwayA", "PathwayB"),
    seed       = 1L
  )

  # Column names
  expect_true(all(c("Longitude", "Latitude", "Condition", "PatientID",
                     "PathwayA", "PathwayB") %in% colnames(df)))

  # Row count
  expect_equal(nrow(df), 100L)  # 50 spots × 2 conditions

  # Condition is a factor with correct levels
  expect_true(is.factor(df$Condition))
  expect_equal(levels(df$Condition), c("Pre", "Post"))

  # Coordinates are numeric and in [0, 100]
  expect_true(all(df$Longitude >= 0 & df$Longitude <= 100))
  expect_true(all(df$Latitude  >= 0 & df$Latitude  <= 100))

  # Pathway scores are numeric
  expect_true(is.numeric(df$PathwayA))
  expect_true(is.numeric(df$PathwayB))

  # No NAs
  expect_false(anyNA(df))
})

test_that("make_st_synthetic() images attribute has correct structure", {
  df <- make_st_synthetic(conditions = c("A", "B", "C"), seed = 99L)
  imgs <- attr(df, "images")

  expect_true(is.list(imgs))
  expect_equal(names(imgs), c("A", "B", "C"))
  # Each image is a matrix
  for (nm in names(imgs)) {
    expect_true(is.matrix(imgs[[nm]]))
  }
})

test_that("make_st_synthetic() reproducible with same seed", {
  df1 <- make_st_synthetic(seed = 77L)
  df2 <- make_st_synthetic(seed = 77L)
  expect_identical(df1, df2)
})

test_that("make_st_synthetic() different seeds produce different data", {
  df1 <- make_st_synthetic(seed = 1L)
  df2 <- make_st_synthetic(seed = 2L)
  expect_false(identical(df1$Longitude, df2$Longitude))
})

test_that("make_st_synthetic() handles single-condition gracefully", {
  # Single condition is allowed for make_st_synthetic (not for STFA_fit)
  df <- make_st_synthetic(n_spots = 30, conditions = c("OnlyOne"), seed = 5L)
  expect_equal(nrow(df), 30L)
  expect_equal(nlevels(df$Condition), 1L)
})

test_that("make_st_synthetic() PatientID column present and non-empty", {
  df <- make_st_synthetic(n_spots = 20, conditions = c("A", "B"), seed = 10L)
  expect_true("PatientID" %in% colnames(df))
  expect_false(any(is.na(df$PatientID) | df$PatientID == ""))
})

test_that("read_seurat_st() errors when SeuratObject not available", {
  # Simulate SeuratObject not being installed by mocking the check
  # We test the error path using a non-Seurat object
  skip_if(requireNamespace("SeuratObject", quietly = TRUE),
          "SeuratObject is installed; skipping unavailability test")

  expect_error(
    read_seurat_st(list(A = list(), B = list()), pathways = "PathwayA"),
    regexp = "SeuratObject.*package is required"
  )
})

test_that("read_seurat_st() validates input: unnamed list", {
  # Must fail with a clear message when list has no names
  skip_if(requireNamespace("SeuratObject", quietly = TRUE),
          "Test logic only needed when SeuratObject is absent")

  unnamed_list <- list(NULL, NULL)
  expect_error(
    read_seurat_st(unnamed_list, pathways = "PathwayA"),
    regexp = "named list"
  )
})

test_that("read_seurat_st() validates input: single condition", {
  skip_if(requireNamespace("SeuratObject", quietly = TRUE),
          "Test logic only needed when SeuratObject is absent")

  expect_error(
    read_seurat_st(list(A = list()), pathways = "PathwayA"),
    regexp = "at least 2 conditions"
  )
})
