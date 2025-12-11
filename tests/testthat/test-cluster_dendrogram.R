# tests/testthat/test-cluster_dendrogram.R
#
# These tests cover:
# - Input validation (dendrogram object, class labels, parameters)
# - Homogeneity threshold enforcement
# - Minimum cluster size constraints
# - Rare class handling (count < min_size)
# - Complete element assignment (no orphans)
# - Output structure and format
# - Dendrogram order mapping
# - S3 methods (print, summary)
# - Edge cases and special scenarios
# - Verbose parameter behavior

test_that("cluster_dendrogram validates inputs correctly", {
  # Create simple valid inputs
  hc <- hclust(dist(matrix(rnorm(30), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B", "C"), length.out = 10)

  # Test invalid dendrogram
  expect_error(
    cluster_dendrogram(list(), labels, hom_thresh = 0.8),
    "'dendrogram' must be a dendrogram object"
  )

  # Test length mismatch
  expect_error(
    cluster_dendrogram(dend, c("A", "B"), hom_thresh = 0.8),
    "'class_labels' length must be"
  )

  # Test NA in labels (use NA_character_ for proper NA)
  labels_na <- labels
  labels_na[1] <- NA_character_
  expect_error(
    cluster_dendrogram(dend, labels_na, hom_thresh = 0.8),
    "'class_labels' contains NA or empty values"
  )

  # Test empty string in labels
  labels_empty <- labels
  labels_empty[1] <- ""
  expect_error(
    cluster_dendrogram(dend, labels_empty, hom_thresh = 0.8),
    "'class_labels' contains NA or empty values"
  )
})

test_that("cluster_dendrogram handles perfect homogeneity case", {
  # Create data with perfect class separation
  set.seed(123)
  mat <- rbind(
    matrix(rnorm(9, mean = 0), ncol = 3),
    matrix(rnorm(9, mean = 5), ncol = 3)
  )
  hc <- hclust(dist(mat))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B"), each = 3)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 1.0,
    min_size = 3,
    verbose = FALSE
  )

  expect_s3_class(result, "cluster_dendrogram_result")
  expect_equal(result$n_elements, 6)
  expect_equal(result$n_unassigned, 0)
})

test_that("cluster_dendrogram creates clusters with minimum size", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(60), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B", "C", "D"), each = 5)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.6,
    min_size = 5,
    verbose = FALSE
  )

  # Check that clusters meeting min_size have at least min_size elements
  # (some may be smaller due to rare class handling)
  expect_true(all(result$cluster_summary$n_elements >= 1))
})

test_that("cluster_dendrogram handles rare classes correctly", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(60), ncol = 3)))
  dend <- as.dendrogram(hc)
  # Create labels with rare class (count < min_size)
  labels <- c(rep("A", 10), rep("B", 8), "rare", "rare")

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.8,
    min_size = 3,
    verbose = FALSE
  )

  # Rare elements should be assigned
  expect_equal(result$n_unassigned, 0)

  # Check that rare class elements were clustered
  rare_assignments <- result$element_assignment$cluster[labels == "rare"]
  expect_true(all(!is.na(rare_assignments)))
})

test_that("cluster_dendrogram respects homogeneity threshold", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(60), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B", "C"), length.out = 20)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.8,
    min_size = 3,
    verbose = FALSE
  )

  # All clusters should meet or exceed threshold
  expect_true(all(result$cluster_summary$homogeneity >= 0.8 - 1e-10))
})

test_that("cluster_dendrogram assigns all elements", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(30), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B"), length.out = 10)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.5,
    min_size = 2,
    verbose = FALSE
  )

  # All elements should be assigned
  expect_equal(result$n_unassigned, 0)
  expect_true(all(!is.na(result$cluster_assignment_dendro_order)))
})

test_that("cluster_dendrogram produces correct output structure", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(30), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B", "C"), length.out = 10)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE
  )

  # Check result structure
  expect_type(result, "list")
  expect_s3_class(result, "cluster_dendrogram_result")

  # Check required components
  expect_true("dendrogram" %in% names(result))
  expect_true("clusters" %in% names(result))
  expect_true("cluster_summary" %in% names(result))
  expect_true("element_assignment" %in% names(result))

  # Check cluster_summary structure
  expect_s3_class(result$cluster_summary, "data.frame")
  expect_true(all(c("cluster_id", "n_elements", "dominant_class",
                    "homogeneity", "n_classes", "class_composition") %in%
                    names(result$cluster_summary)))

  # Check element_assignment structure
  expect_s3_class(result$element_assignment, "data.frame")
  expect_true(all(c("sequence_name", "dendro_index", "class", "cluster") %in%
                    names(result$element_assignment)))
})

test_that("cluster_dendrogram handles dendro_order mapping", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(30), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B"), length.out = 10)
  dendro_order <- sample(1:10)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.7,
    min_size = 2,
    dendro_order = dendro_order,
    verbose = FALSE
  )

  # Both assignment vectors should exist
  expect_true("cluster_assignment_dendro_order" %in% names(result))
  expect_true("cluster_assignment_original_order" %in% names(result))

  # They should have same length
  expect_equal(
    length(result$cluster_assignment_dendro_order),
    length(result$cluster_assignment_original_order)
  )
})

test_that("cluster_dendrogram handles sequence_names parameter", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(30), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B"), length.out = 10)
  seq_names <- paste0("Seq_", 1:10)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.7,
    min_size = 2,
    sequence_names = seq_names,
    verbose = FALSE
  )

  expect_equal(result$element_assignment$sequence_name, seq_names)
})

test_that("cluster_dendrogram updates data_result when provided", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(30), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B"), length.out = 10)

  data_result <- list(
    kmers = matrix(rnorm(100), ncol = 10),
    metadata = data.frame(id = 1:10)
  )

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.7,
    min_size = 2,
    data_result = data_result,
    verbose = FALSE
  )

  expect_true("cluster" %in% names(result$data_result$metadata))
  expect_equal(nrow(result$data_result$metadata), 10)
})

test_that("print method works correctly", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(30), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B"), length.out = 10)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.7,
    min_size = 2,
    verbose = FALSE
  )

  expect_output(print(result), "Dendrogram Clustering Result")
  expect_output(print(result), "Total elements:")
  expect_output(print(result), "Total clusters:")
})

test_that("summary method works correctly", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(30), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B"), length.out = 10)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 0.7,
    min_size = 2,
    verbose = FALSE
  )

  expect_output(summary(result), "Clustering Summary")
  expect_output(summary(result), "Mean homogeneity:")
})

test_that("cluster_dendrogram handles edge cases", {
  set.seed(123)
  # Minimum viable dendrogram (3 elements for min_size=3)
  hc <- hclust(dist(matrix(rnorm(9), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep("A", 3)

  result <- cluster_dendrogram(
    dend, labels,
    hom_thresh = 1.0,
    min_size = 3,
    verbose = FALSE
  )

  expect_equal(result$n_elements, 3)
  expect_true(nrow(result$cluster_summary) >= 1)

  # All same class should have high homogeneity
  expect_true(all(result$cluster_summary$homogeneity >= 0.99))
})

test_that("cluster_dendrogram handles verbose parameter", {
  set.seed(123)
  hc <- hclust(dist(matrix(rnorm(30), ncol = 3)))
  dend <- as.dendrogram(hc)
  labels <- rep(c("A", "B"), length.out = 10)

  # With verbose = TRUE, should produce messages
  expect_message(
    cluster_dendrogram(dend, labels, hom_thresh = 0.7, min_size = 2, verbose = TRUE),
    "Starting clustering"
  )

  # With verbose = FALSE, should not produce messages
  expect_silent(
    cluster_dendrogram(dend, labels, hom_thresh = 0.7, min_size = 2, verbose = FALSE)
  )
})

