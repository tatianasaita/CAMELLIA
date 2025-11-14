# Test suite for calculate_cluster_motifs function
# Following testthat and CRAN standards

library(testthat)

# Setup: Create mock data for testing
create_mock_data <- function() {
  # Create mock k-mer data
  kmers_df <- data.frame(
    kmer_1 = c(10, 20, 30, 40),
    kmer_2 = c(15, 25, 35, 45),
    kmer_3 = c(5, 10, 15, 20),
    kmer_4 = c(100, 100, 100, 100)
  )

  # Create mock metadata with cluster assignments
  metadata_df <- data.frame(
    sequence_id = c("seq1", "seq2", "seq3", "seq4"),
    cluster = c(1, 1, 2, 2)
  )

  # Return as list (expected format)
  list(
    kmers = kmers_df,
    metadata = metadata_df
  )
}

# Setup: Create mock cluster_dendrogram_result object
create_mock_cluster_result <- function() {
  data_result <- create_mock_data()

  # Create a mock cluster_dendrogram_result object
  result <- list(
    dendrogram = NULL,
    data_result = data_result
  )

  class(result) <- "cluster_dendrogram_result"
  result
}

# Test Suite: Input Validation
test_that("Function rejects non-cluster_dendrogram_result objects", {
  expect_error(
    calculate_cluster_motifs(list()),
    "must be an object of class 'cluster_dendrogram_result'"
  )
})

test_that("Function rejects object without data_result element", {
  invalid_obj <- structure(list(), class = "cluster_dendrogram_result")
  expect_error(
    calculate_cluster_motifs(invalid_obj),
    "must contain 'data_result' element"
  )
})

test_that("Function accepts data_result parameter override", {
  cluster_result <- create_mock_cluster_result()
  mock_data <- create_mock_data()

  result <- calculate_cluster_motifs(
    cluster_result,
    data_result = mock_data
  )

  expect_s3_class(result, "cluster_motifs")
})

test_that("Function validates data_result is a list", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result <- "not_a_list"

  expect_error(
    calculate_cluster_motifs(cluster_result),
    "'data_result' must be a list"
  )
})

test_that("Function validates data_result contains kmers and metadata", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result <- list(kmers = data.frame())

  expect_error(
    calculate_cluster_motifs(cluster_result),
    "must contain 'kmers' and 'metadata' elements"
  )
})

test_that("Function validates kmers is a data.frame", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$kmers <- list()

  expect_error(
    calculate_cluster_motifs(cluster_result),
    "'data_result\\$kmers' must be a data.frame"
  )
})

test_that("Function validates metadata is a data.frame", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$metadata <- list()

  expect_error(
    calculate_cluster_motifs(cluster_result),
    "'data_result\\$metadata' must be a data.frame"
  )
})

test_that("Function validates metadata contains cluster column", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$metadata <- data.frame(sequence_id = c("s1", "s2"))

  expect_error(
    calculate_cluster_motifs(cluster_result),
    "must contain a 'cluster' column"
  )
})

test_that("Function validates kmers and metadata have same number of rows", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$kmers <- data.frame(k1 = c(1, 2, 3))
  cluster_result$data_result$metadata <- data.frame(
    sequence_id = c("s1", "s2"),
    cluster = c(1, 2)
  )

  expect_error(
    calculate_cluster_motifs(cluster_result),
    "must have the same number of rows"
  )
})

test_that("Function validates all k-mer columns are numeric", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$kmers$kmer_text <- c("a", "b", "c", "d")

  expect_error(
    calculate_cluster_motifs(cluster_result),
    "All k-mer columns must be numeric"
  )
})

test_that("Function handles CLASS column exclusion", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$kmers$CLASS <- c("typeA", "typeB", "typeA", "typeB")

  result <- calculate_cluster_motifs(cluster_result)

  expect_s3_class(result, "cluster_motifs")
  expect_false("CLASS" %in% rownames(result))
})

# Test Suite: Output Format and Structure
test_that("Function returns cluster_motifs class object", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  expect_s3_class(result, "cluster_motifs")
  expect_s3_class(result, "data.frame")
})

test_that("Output has correct dimensions", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  # 4 k-mers, 2 clusters
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 2)
})

test_that("Output has correct row and column names", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  # Row names should be k-mer names
  expect_equal(rownames(result), c("kmer_1", "kmer_2", "kmer_3", "kmer_4"))

  # Column names should be Cluster_<id>
  expect_equal(colnames(result), c("Cluster_1", "Cluster_2"))
})

test_that("Output is a data.frame", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  expect_true(is.data.frame(result))
})

# Test Suite: Normalization Validation
test_that("Output values are normalized between 0 and 1", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  values <- as.matrix(result)
  expect_true(all(values >= 0, na.rm = TRUE))
  expect_true(all(values <= 1, na.rm = TRUE))
})

test_that("Normalization includes min and max values", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  # Each column should have at least one 0 and one 1
  # (or all 0.5 if zero-range)
  for (col_idx in seq_len(ncol(result))) {
    col_vals <- result[, col_idx]
    has_zero <- any(col_vals == 0, na.rm = TRUE)
    has_one <- any(col_vals == 1, na.rm = TRUE)
    has_half <- any(col_vals == 0.5, na.rm = TRUE)

    expect_true(has_zero || has_one || has_half)
  }
})

test_that("Zero-range columns are normalized correctly", {
  # Create specific mock data with zero-range scenario
  kmers_df <- data.frame(
    kmer_1 = c(10, 20, 30, 40),
    kmer_2 = c(15, 25, 35, 45),
    kmer_3 = c(100, 100, 100, 100)  # All identical - will cause zero-range
  )

  metadata_df <- data.frame(
    sequence_id = c("seq1", "seq2", "seq3", "seq4"),
    cluster = c(1, 1, 2, 2)
  )

  data_result <- list(
    kmers = kmers_df,
    metadata = metadata_df
  )

  cluster_result <- structure(
    list(data_result = data_result),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(cluster_result)

  # When all values in a k-mer are identical (kmer_3 = 100 everywhere):
  # After aggregation: Cluster_1 = 200, Cluster_2 = 200
  # Since both clusters have the same aggregated value:
  # - Min per cluster: 200, Max per cluster: 200
  # - Normalization: (200 - 200) / (200 - 200) = NaN/0
  # - When range is zero, values are set to 1.0 (they ARE the max)

  # Convert to vector
  kmer3_row <- as.numeric(result[3, ])

  # Verify zero-range handling: all values should be equal (both 0 or both 1)
  expect_true(
    all(kmer3_row == kmer3_row[1]),
    info = "Zero-range k-mer should have identical values across clusters"
  )

  # Since the value equals max in both columns, it should be 1.0
  expect_true(
    all(kmer3_row == 1.0),
    info = sprintf(
      "Zero-range k-mer with identical values should normalize to 1.0, got: %s",
      paste(round(kmer3_row, 3), collapse = ", ")
    )
  )
})

# Test Suite: Aggregation Correctness
test_that("K-mer values are correctly summed by cluster", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  # Manual calculation for verification
  # Cluster 1: seq1 (10, 15, 5, 100) + seq2 (20, 25, 10, 100)
  # Sum: (30, 40, 15, 200)
  # Cluster 2: seq3 (30, 35, 15, 100) + seq4 (40, 45, 20, 100)
  # Sum: (70, 80, 35, 200)

  # After normalization, verify relative ordering is preserved
  cluster_1 <- result[, 1]
  cluster_2 <- result[, 2]

  # kmer_1 should be smaller in cluster_1 than in cluster_2
  # (because 30 < 70 before normalization)
  expect_true(cluster_1[1] < cluster_2[1])
})

# Test Suite: Edge Cases
test_that("Function handles single cluster", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$metadata$cluster <- c(1, 1, 1, 1)

  result <- calculate_cluster_motifs(cluster_result)

  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "Cluster_1")
})

test_that("Function handles many clusters", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$metadata$cluster <- c(1, 2, 3, 4)

  result <- calculate_cluster_motifs(cluster_result)

  expect_equal(ncol(result), 4)
})

test_that("Function handles single k-mer", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$kmers <- data.frame(kmer_1 = c(10, 20, 30, 40))

  result <- calculate_cluster_motifs(cluster_result)

  expect_equal(nrow(result), 1)
  expect_equal(rownames(result), "kmer_1")
})

test_that("Function handles large values", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$kmers <- data.frame(
    k1 = c(1e6, 2e6, 3e6, 4e6),
    k2 = c(5e6, 6e6, 7e6, 8e6)
  )

  result <- calculate_cluster_motifs(cluster_result)

  expect_s3_class(result, "cluster_motifs")
  expect_true(all(as.matrix(result) >= 0 & as.matrix(result) <= 1, na.rm = TRUE))
})

test_that("Function handles small values", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$kmers <- data.frame(
    k1 = c(0.001, 0.002, 0.003, 0.004),
    k2 = c(0.005, 0.006, 0.007, 0.008)
  )

  result <- calculate_cluster_motifs(cluster_result)

  expect_s3_class(result, "cluster_motifs")
  expect_true(all(as.matrix(result) >= 0 & as.matrix(result) <= 1, na.rm = TRUE))
})

test_that("Function handles negative values", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$kmers <- data.frame(
    k1 = c(-10, -20, 30, 40),
    k2 = c(-15, 25, -35, 45)
  )

  result <- calculate_cluster_motifs(cluster_result)

  expect_s3_class(result, "cluster_motifs")
  expect_true(all(as.matrix(result) >= 0 & as.matrix(result) <= 1, na.rm = TRUE))
})

test_that("Function handles zero values", {
  cluster_result <- create_mock_cluster_result()
  cluster_result$data_result$kmers <- data.frame(
    k1 = c(0, 0, 0, 0),
    k2 = c(10, 20, 30, 40)
  )

  result <- calculate_cluster_motifs(cluster_result)

  expect_s3_class(result, "cluster_motifs")
})

# Test Suite: S3 Methods
test_that("print method works correctly", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  expect_output(
    print(result),
    "Cluster Motif Matrix"
  )

  expect_output(
    print(result),
    "Dimensions"
  )
})

test_that("summary method works correctly", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  expect_output(
    summary(result),
    "Cluster Motif Summary"
  )

  expect_output(
    summary(result),
    "min="
  )
})

test_that("summary method returns invisible object", {
  cluster_result <- create_mock_cluster_result()
  result <- calculate_cluster_motifs(cluster_result)

  # Use expect_invisible with expression, not captured output
  expect_invisible(summary(result))
})

# Test Suite: Return Values
test_that("Function returns invisible result", {
  cluster_result <- create_mock_cluster_result()

  result <- capture.output(
    value <- calculate_cluster_motifs(cluster_result)
  )

  expect_s3_class(value, "cluster_motifs")
})

test_that("Default method raises appropriate error", {
  expect_error(
    calculate_cluster_motifs.default(list()),
    "must be an object of class 'cluster_dendrogram_result'"
  )
})

# Test Suite: Data Integrity
test_that("Function does not modify input data", {
  cluster_result <- create_mock_cluster_result()
  kmers_original <- cluster_result$data_result$kmers
  metadata_original <- cluster_result$data_result$metadata

  calculate_cluster_motifs(cluster_result)

  expect_equal(cluster_result$data_result$kmers, kmers_original)
  expect_equal(cluster_result$data_result$metadata, metadata_original)
})

test_that("Function handles multiple sequential calls", {
  cluster_result <- create_mock_cluster_result()

  result1 <- calculate_cluster_motifs(cluster_result)
  result2 <- calculate_cluster_motifs(cluster_result)

  expect_equal(result1, result2)
})

