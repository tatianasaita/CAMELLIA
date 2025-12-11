# tests/testthat/test-kmer_analysis.R
#
# Essential tests for kmer_analysis function following CRAN standards
# These tests cover:
# - Input validation (cluster_result, motif_matrix, threshold)
# - Core functionality (unique motifs identification, frequency ranking)
# - Validation data handling (external validation datasets)
# - Threshold behavior (motif presence/absence logic)
# - S3 methods (print, summary)
# - Edge cases (empty results, single cluster/class scenarios)
# - Data structure integrity (output components, mappings)

library(testthat)

test_that("kmer_analysis validates inputs", {
  # Invalid cluster_result
  expect_error(
    kmer_analysis(cluster_result = list()),
    "cluster_summary"
  )

  # Invalid threshold
  expect_error(
    kmer_analysis(
      cluster_result = list(cluster_summary = data.frame()),
      threshold = "invalid"
    )
  )
})

test_that("kmer_analysis works with minimal cluster_result", {
  # Create minimal valid cluster_result
  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = 1:3,
      dominant_class = c("A", "B", "A"),
      stringsAsFactors = FALSE
    ),
    data_result = NULL
  )

  # Create motif matrix
  motif_matrix <- data.frame(
    Cluster_1 = c(0.8, 0.1, 0.0),
    Cluster_2 = c(0.0, 0.9, 0.1),
    Cluster_3 = c(0.7, 0.0, 0.0),
    row.names = c("AAA", "TTT", "GGG")
  )

  result <- kmer_analysis(cluster_result, motif_matrix, threshold = 0)

  expect_s3_class(result, "kmer_analysis_result")
  expect_true("unique_cluster_motifs" %in% names(result))
  expect_true("cluster_frequency_ranking" %in% names(result))
  expect_true("unique_class_motifs" %in% names(result))
  expect_true("class_frequency_matrix" %in% names(result))
  expect_false(result$has_external_validation)
})

test_that("kmer_analysis identifies unique cluster motifs", {
  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = 1:2,
      dominant_class = c("A", "B"),
      stringsAsFactors = FALSE
    )
  )

  motif_matrix <- data.frame(
    Cluster_1 = c(0.9, 0.0),
    Cluster_2 = c(0.0, 0.8),
    row.names = c("AAA", "TTT")
  )

  result <- kmer_analysis(cluster_result, motif_matrix, threshold = 0)

  expect_equal(nrow(result$unique_cluster_motifs), 2)
  expect_true(all(c("motif", "cluster", "class", "value") %in%
                    colnames(result$unique_cluster_motifs)))
})

test_that("kmer_analysis handles validation data", {
  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = 1:2,
      dominant_class = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    classification_result = list(
      validation_data = data.frame(
        CLASS = c("A", "B", "A"),
        AAA = c(1, 0, 1),
        TTT = c(0, 1, 0),
        stringsAsFactors = FALSE
      )
    )
  )

  motif_matrix <- data.frame(
    Cluster_1 = c(0.8, 0.1),
    Cluster_2 = c(0.1, 0.9),
    row.names = c("AAA", "TTT")
  )

  result <- kmer_analysis(cluster_result, motif_matrix)

  expect_true(result$has_external_validation)
  expect_false(is.null(result$validation))
  expect_true("unique_class_motifs" %in% names(result$validation))
  expect_true("class_frequency_ranking" %in% names(result$validation))
  expect_true("motifs_by_class_rank" %in% names(result$validation))
})

test_that("kmer_analysis threshold works correctly", {
  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = 1:2,
      dominant_class = c("A", "B"),
      stringsAsFactors = FALSE
    )
  )

  motif_matrix <- data.frame(
    Cluster_1 = c(0.5, 0.1),
    Cluster_2 = c(0.1, 0.6),
    row.names = c("AAA", "TTT")
  )

  result_low <- kmer_analysis(cluster_result, motif_matrix, threshold = 0)
  result_high <- kmer_analysis(cluster_result, motif_matrix, threshold = 0.3)

  # Com threshold mais alto, mais motifs são considerados únicos
  # porque valores baixos são filtrados
  expect_lte(nrow(result_low$unique_cluster_motifs),
             nrow(result_high$unique_cluster_motifs))
})

test_that("kmer_analysis print method works", {
  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = 1,
      dominant_class = "A",
      stringsAsFactors = FALSE
    )
  )

  motif_matrix <- data.frame(
    Cluster_1 = c(0.8),
    row.names = c("AAA")
  )

  result <- kmer_analysis(cluster_result, motif_matrix)

  expect_output(print(result), "K-mer Analysis Results")
  expect_output(print(result), "TRAINING DATA")
})

test_that("kmer_analysis summary method works", {
  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = 1:2,
      dominant_class = c("A", "B"),
      stringsAsFactors = FALSE
    )
  )

  motif_matrix <- data.frame(
    Cluster_1 = c(0.8, 0.1, 0.2),
    Cluster_2 = c(0.1, 0.9, 0.3),
    row.names = c("AAA", "TTT", "GGG")
  )

  result <- kmer_analysis(cluster_result, motif_matrix)

  expect_output(summary(result), "TRAINING SUMMARY")
  expect_output(summary(result), "most frequent motifs")
})

test_that("kmer_analysis handles empty results gracefully", {
  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = 1,
      dominant_class = "A",
      stringsAsFactors = FALSE
    )
  )

  # All motifs present in all clusters
  motif_matrix <- data.frame(
    Cluster_1 = c(0.8, 0.9),
    row.names = c("AAA", "TTT")
  )

  result <- kmer_analysis(cluster_result, motif_matrix, threshold = 0.5)

  expect_s3_class(result, "kmer_analysis_result")
  expect_true(is.data.frame(result$unique_cluster_motifs))
  expect_true(is.data.frame(result$unique_class_motifs))
})

test_that("kmer_analysis cluster_to_class mapping is correct", {
  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = c(1, 2, 3),
      dominant_class = c("A", "B", "A"),
      stringsAsFactors = FALSE
    )
  )

  motif_matrix <- data.frame(
    Cluster_1 = c(0.8),
    Cluster_2 = c(0.7),
    Cluster_3 = c(0.6),
    row.names = c("AAA")
  )

  result <- kmer_analysis(cluster_result, motif_matrix)

  expected_mapping <- c(Cluster_1 = "A", Cluster_2 = "B", Cluster_3 = "A")

  expect_equal(result$cluster_to_class, expected_mapping)
})

