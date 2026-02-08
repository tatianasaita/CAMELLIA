# tests/testthat/test-select_motifs.R

library('testthat')

test_that("select_motifs returns correct structure for valid inputs", {
  # Setup test data
  motif_cluster <- data.frame(
    Cluster_1 = c(0.8, 0.6, 0.4, 0.2),
    Cluster_2 = c(0.3, 0.7, 0.5, 0.1),
    Cluster_3 = c(0.2, 0.4, 0.9, 0.6),
    row.names = c("AAAA", "AAAT", "AATG", "AACG")
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:3,
    dominant_class = c("ClassA", "ClassA", "ClassB")
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  # Test Case 2: k <= n < m
  result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)
  
  expect_s3_class(result, "select_motifs")
  expect_type(result, "list")
  expect_equal(attr(result, "n"), 2L)
  expect_equal(attr(result, "k"), 2L)
  expect_equal(attr(result, "m"), 3L)
  expect_equal(attr(result, "case"), "CASE_2")
})

test_that("select_motifs handles Case 1: n < k", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.8, 0.6),
    Cluster_2 = c(0.3, 0.7),
    row.names = c("AAAA", "AAAT")
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:2,
    dominant_class = c("ClassA", "ClassB")
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  result <- select_motifs(motif_cluster, cluster_result, n = 1, verbose = FALSE)
  
  expect_s3_class(result, "empty_result")
  expect_equal(attr(result, "case"), "CASE_1")
  expect_equal(attr(result, "reason"), "n < k: insufficient motifs for all classes")
  expect_length(result, 0)
})

test_that("select_motifs handles Case 3: n >= m", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.8, 0.6, 0.4, 0.2, 0.1),
    Cluster_2 = c(0.3, 0.7, 0.5, 0.1, 0.2),
    row.names = c("AAAA", "AAAT", "AATG", "AACG", "ACGT")
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:2,
    dominant_class = c("ClassA", "ClassB")
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  result <- select_motifs(motif_cluster, cluster_result, n = 4, verbose = FALSE)
  
  expect_s3_class(result, "select_motifs")
  expect_equal(attr(result, "case"), "CASE_3")
  expect_equal(sum(lengths(result)), 4)
})

test_that("select_motifs validates input types", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.8, 0.6),
    row.names = c("AAAA", "AAAT")
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1,
    dominant_class = "ClassA"
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  # Test integer conversion
  result <- select_motifs(motif_cluster, cluster_result, n = 1.5, verbose = FALSE)
  expect_equal(attr(result, "n"), 1L)
})

test_that("select_motifs ensures no duplicate motifs", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.8, 0.7),
    Cluster_2 = c(0.6, 0.5, 0.4),
    row.names = c("AAAA", "AAAT", "AATG")
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:2,
    dominant_class = c("ClassA", "ClassB")
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)
  
  all_motifs <- unlist(result, use.names = FALSE)
  expect_equal(length(all_motifs), length(unique(all_motifs)))
})

test_that("select_motifs returns correct attributes", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.8, 0.6, 0.4),
    Cluster_2 = c(0.3, 0.7, 0.5),
    row.names = c("AAAA", "AAAT", "AATG")
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:2,
    dominant_class = c("ClassA", "ClassB")
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)
  
  expect_true("classes_order" %in% names(attributes(result)))
  expect_true("motifs_per_class" %in% names(attributes(result)))
  expect_type(attr(result, "classes_order"), "character")
  expect_type(attr(result, "motifs_per_class"), "integer")
})

