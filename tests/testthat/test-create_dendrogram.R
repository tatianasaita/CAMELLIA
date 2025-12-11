# ==============================================================================
# File: tests/testthat/test_create_dendrogram.R
# ==============================================================================

skip_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    skip(paste("Package", pkg, "not available"))
  }
}

test_that("create_dendrogram requires CLASS column", {
  data <- data.frame(
    x = rnorm(10),
    y = rnorm(10)
  )

  expect_error(
    create_dendrogram(data),
    "'data' must contain a 'CLASS' column"
  )
})

test_that("create_dendrogram requires at least 2 rows", {
  data <- data.frame(
    x = 1,
    CLASS = "A"
  )

  expect_error(
    create_dendrogram(data),
    "'data' must contain at least 2 rows"
  )
})

test_that("create_dendrogram requires numeric columns", {
  data <- data.frame(
    x = c("a", "b"),
    CLASS = c("A", "B")
  )

  expect_error(
    create_dendrogram(data),
    "'data' must contain numeric columns"
  )
})

test_that("create_dendrogram validates sequence_names length", {
  data <- data.frame(
    x = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  expect_error(
    create_dendrogram(data, sequence_names = c("seq1", "seq2")),
    "'sequence_names' must have length 10"
  )
})

test_that("create_dendrogram rejects NA values", {
  data <- data.frame(
    x = c(rnorm(9), NA),
    CLASS = rep(c("A", "B"), 5)
  )

  expect_error(
    create_dendrogram(data),
    "contains NA values"
  )
})

test_that("create_dendrogram rejects Inf values", {
  data <- data.frame(
    x = c(rnorm(9), Inf),
    CLASS = rep(c("A", "B"), 5)
  )

  expect_error(
    create_dendrogram(data),
    "contains Inf values"
  )
})

test_that("create_dendrogram validates dist_method", {
  data <- data.frame(
    x = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  expect_error(
    create_dendrogram(data, dist_method = "invalid_method"),
    "'dist_method' must be one of"
  )
})

test_that("create_dendrogram validates hclust_method", {
  data <- data.frame(
    x = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  expect_error(
    create_dendrogram(data, hclust_method = "invalid_method"),
    "'hclust_method' must be one of"
  )
})

test_that("create_dendrogram validates show_labels", {
  data <- data.frame(
    x = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  expect_error(
    create_dendrogram(data, show_labels = "yes"),
    "'show_labels' must be a logical value"
  )
})

test_that("create_dendrogram validates plot_title", {
  data <- data.frame(
    x = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  expect_error(
    create_dendrogram(data, plot_title = c("Title1", "Title2")),
    "'plot_title' must be a single character string"
  )
})

test_that("create_dendrogram produces correct output structure", {
  set.seed(123)
  data <- data.frame(
    x = rnorm(10),
    y = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  # ← FIXED: Suppress dendextend warning
  result <- suppressWarnings(
    create_dendrogram(data, plot_title = "Test Dendrogram")
  )

  expect_s3_class(result, "dendrogram_result")
  expect_named(result, c("dendrogram", "hclust", "order", "labels",
                         "sequence_names", "colors", "base_colors",
                         "height", "n_samples", "n_classes"))

  expect_s3_class(result$dendrogram, "dendrogram")
  expect_s3_class(result$hclust, "hclust")
  expect_true(is.integer(result$order))
  expect_true(is.character(result$labels))
  expect_true(is.character(result$sequence_names))
  expect_true(is.character(result$colors))
  expect_true(is.character(result$base_colors))
  expect_true(is.numeric(result$height))
  expect_true(is.integer(result$n_samples))
  expect_true(is.integer(result$n_classes))
})

test_that("create_dendrogram colors mapped to classes", {
  set.seed(123)
  data <- data.frame(
    x = rnorm(10),
    y = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  # ← FIXED: Suppress warning
  result <- suppressWarnings(create_dendrogram(data))

  expect_named(result$base_colors, c("A", "B"))
  expect_equal(length(result$colors), 10)

  for (i in seq_along(result$colors)) {
    class_label <- result$labels[i]
    expected_color <- result$base_colors[class_label]
    expect_equal(result$colors[i], expected_color)
  }
})

test_that("create_dendrogram with sequence_names", {
  set.seed(123)
  data <- data.frame(
    x = rnorm(10),
    y = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  seq_names <- paste0("seq_", 1:10)

  # ← FIXED: Suppress warning
  result <- suppressWarnings(
    create_dendrogram(data, sequence_names = seq_names)
  )

  expect_equal(length(result$sequence_names), 10)
  expect_true(all(grepl("^seq_", result$sequence_names)))
})

test_that("create_dendrogram different distance methods", {
  set.seed(123)
  data <- data.frame(
    x = rnorm(10),
    y = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  dist_methods <- c("euclidean", "manhattan", "maximum")

  for (method in dist_methods) {
    # ← FIXED: Suppress warning
    result <- suppressWarnings(
      create_dendrogram(data, dist_method = method)
    )
    expect_s3_class(result, "dendrogram_result")
    expect_true(!is.null(result$dendrogram))
  }
})

test_that("create_dendrogram different clustering methods", {
  set.seed(123)
  data <- data.frame(
    x = rnorm(10),
    y = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  hclust_methods <- c("ward.D2", "single", "complete", "average")

  for (method in hclust_methods) {
    # ← FIXED: Suppress warning
    result <- suppressWarnings(
      create_dendrogram(data, hclust_method = method)
    )
    expect_s3_class(result, "dendrogram_result")
    expect_true(!is.null(result$dendrogram))
  }
})

test_that("create_dendrogram with multiple classes", {
  set.seed(123)
  data <- data.frame(
    x = rnorm(30),
    y = rnorm(30),
    CLASS = rep(c("A", "B", "C"), 10)
  )

  # ← FIXED: Suppress warning
  result <- suppressWarnings(create_dendrogram(data))

  expect_equal(result$n_classes, 3)
  expect_equal(length(result$base_colors), 3)
  expect_named(result$base_colors, c("A", "B", "C"))
})

test_that("print.dendrogram_result displays correctly", {
  set.seed(123)
  data <- data.frame(
    x = rnorm(10),
    y = rnorm(10),
    CLASS = rep(c("A", "B"), 5)
  )

  # ← FIXED: Suppress warning
  result <- suppressWarnings(create_dendrogram(data))

  expect_output(print(result), "Dendrogram Analysis Result")
  expect_output(print(result), "Number of samples")
  expect_output(print(result), "Number of classes")
  expect_output(print(result), "Color mapping by class")
})
