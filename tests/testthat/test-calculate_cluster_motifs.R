# tests/testthat/test-calculate_cluster_motifs.R

library(testthat)

test_that("calculate_cluster_motifs returns correct structure", {
  # Setup mock data
  kmers_df <- data.frame(
    kmer1 = c(10, 20, 30, 40),
    kmer2 = c(15, 25, 35, 45),
    kmer3 = c(5, 10, 15, 20)
  )

  data_result <- list(kmers = kmers_df)

  result_cluster_dendrogram <- list(
    cluster_assignment_dendro_order = c(1, 1, 2, 2),
    data_result = data_result
  )

  # Run function
  result <- calculate_cluster_motifs(result_cluster_dendrogram)

  # Test output structure
  expect_s3_class(result, "cluster_motifs")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)  # 3 k-mers
  expect_equal(ncol(result), 2)  # 2 clusters
  expect_named(result, c("Cluster_1", "Cluster_2"))
})

test_that("calculate_cluster_motifs normalizes correctly", {
  # Setup mock data
  kmers_df <- data.frame(
    kmer1 = c(10, 20, 30, 40),
    kmer2 = c(0, 0, 100, 100)
  )

  data_result <- list(kmers = kmers_df)

  result_cluster_dendrogram <- list(
    cluster_assignment_dendro_order = c(1, 1, 2, 2),
    data_result = data_result
  )

  # Run function
  result <- calculate_cluster_motifs(result_cluster_dendrogram)

  # Test normalization (values should be between 0 and 1)
  expect_true(all(result >= 0 & result <= 1))

  # Each column should have at least one 0 and one 1 (min-max normalization)
  expect_true(all(apply(result, 2, min) == 0))
  expect_true(all(apply(result, 2, max) == 1))
})

test_that("calculate_cluster_motifs handles CLASS column correctly", {
  # Setup mock data with CLASS column
  kmers_df <- data.frame(
    kmer1 = c(10, 20, 30),
    kmer2 = c(15, 25, 35),
    CLASS = c("A", "B", "A")
  )

  data_result <- list(kmers = kmers_df)

  result_cluster_dendrogram <- list(
    cluster_assignment_dendro_order = c(1, 2, 1),
    data_result = data_result
  )

  # Run function
  result <- calculate_cluster_motifs(result_cluster_dendrogram)

  # Test that CLASS column is excluded
  expect_equal(nrow(result), 2)  # Only 2 k-mers (CLASS excluded)
  expect_equal(rownames(result), c("kmer1", "kmer2"))
})

test_that("calculate_cluster_motifs handles single cluster", {
  # Setup mock data
  kmers_df <- data.frame(
    kmer1 = c(10, 20, 30),
    kmer2 = c(5, 10, 15)
  )

  data_result <- list(kmers = kmers_df)

  result_cluster_dendrogram <- list(
    cluster_assignment_dendro_order = c(1, 1, 1),
    data_result = data_result
  )

  # Run function
  result <- calculate_cluster_motifs(result_cluster_dendrogram)

  # Test single cluster output
  expect_equal(ncol(result), 1)
  expect_named(result, "Cluster_1")
  expect_equal(nrow(result), 2)
})

test_that("calculate_cluster_motifs row names are k-mer names", {
  # Setup mock data
  kmers_df <- data.frame(
    AAA = c(10, 20),
    CCC = c(15, 25),
    GGG = c(5, 10)
  )

  data_result <- list(kmers = kmers_df)

  result_cluster_dendrogram <- list(
    cluster_assignment_dendro_order = c(1, 2),
    data_result = data_result
  )

  # Run function
  result <- calculate_cluster_motifs(result_cluster_dendrogram)

  # Test row names
  expect_equal(rownames(result), c("AAA", "CCC", "GGG"))
})

