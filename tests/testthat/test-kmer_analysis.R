# tests/testthat/test-kmer_analysis.R

library(testthat)

test_that("kmer_analysis returns correct structure", {
  skip_on_cran()
  
  # Create minimal mock data
  motif_matrix <- matrix(
    c(10, 0, 5,
      0, 8, 0,
      3, 3, 3),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      c("motif_A", "motif_B", "motif_C"),
      c("Cluster_1", "Cluster_2", "Cluster_3")
    )
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:3,
    dominant_class = c("Class_A", "Class_B", "Class_A"),
    stringsAsFactors = FALSE
  )
  
  cluster_result <- list(
    cluster_summary = cluster_summary,
    data_result = NULL,
    classification_result = NULL
  )
  
  # Run function
  result <- kmer_analysis(
    cluster_result = cluster_result,
    motif_matrix = motif_matrix,
    threshold = 0
  )
  
  # Test class
  expect_s3_class(result, "kmer_analysis_result")
  expect_type(result, "list")
  
  # Test essential components
  expect_named(result, c(
    "unique_cluster_motifs",
    "cluster_frequency_ranking",
    "unique_class_motifs",
    "class_frequency_matrix",
    "motifs_by_class_rank",
    "cluster_to_class",
    "has_external_validation",
    "validation"
  ))
  
  # Test data frames structure
  expect_s3_class(result$unique_cluster_motifs, "data.frame")
  expect_s3_class(result$cluster_frequency_ranking, "data.frame")
  expect_s3_class(result$unique_class_motifs, "data.frame")
  expect_s3_class(result$class_frequency_matrix, "data.frame")
  
  # Test validation flag
  expect_type(result$has_external_validation, "logical")
  expect_false(result$has_external_validation)
  expect_null(result$validation)
})

test_that("kmer_analysis identifies unique cluster motifs correctly", {
  skip_on_cran()
  
  # Matrix with clear unique motifs
  motif_matrix <- matrix(
    c(10, 0, 0,   # motif_A unique to Cluster_1
      0, 15, 0,   # motif_B unique to Cluster_2
      5, 5, 5),   # motif_C in all clusters
    nrow = 3, byrow = TRUE,
    dimnames = list(
      c("motif_A", "motif_B", "motif_C"),
      c("Cluster_1", "Cluster_2", "Cluster_3")
    )
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:3,
    dominant_class = c("Class_A", "Class_B", "Class_A"),
    stringsAsFactors = FALSE
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  result <- kmer_analysis(
    cluster_result = cluster_result,
    motif_matrix = motif_matrix,
    threshold = 0
  )
  
  # Should identify 2 unique motifs (A and B)
  expect_equal(nrow(result$unique_cluster_motifs), 2)
  expect_true("motif_A" %in% result$unique_cluster_motifs$motif)
  expect_true("motif_B" %in% result$unique_cluster_motifs$motif)
  expect_false("motif_C" %in% result$unique_cluster_motifs$motif)
  
  # Check columns
  expect_named(result$unique_cluster_motifs, 
               c("motif", "cluster", "class", "value"))
})

test_that("kmer_analysis identifies unique class motifs correctly", {
  skip_on_cran()
  
  motif_matrix <- matrix(
    c(10, 0, 0, 0,   # motif_A only in Class_A clusters
      0, 0, 8, 7,    # motif_B only in Class_B clusters
      5, 5, 5, 5),   # motif_C in both classes
    nrow = 3, byrow = TRUE,
    dimnames = list(
      c("motif_A", "motif_B", "motif_C"),
      c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4")
    )
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:4,
    dominant_class = c("Class_A", "Class_A", "Class_B", "Class_B"),
    stringsAsFactors = FALSE
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  result <- kmer_analysis(
    cluster_result = cluster_result,
    motif_matrix = motif_matrix,
    threshold = 0
  )
  
  # Should identify 2 unique class motifs (A and B)
  expect_equal(nrow(result$unique_class_motifs), 2)
  expect_true("motif_A" %in% result$unique_class_motifs$motif)
  expect_true("motif_B" %in% result$unique_class_motifs$motif)
  
  # Check columns
  expect_named(result$unique_class_motifs,
               c("motif", "class", "n_clusters_with_motif", "total_class_clusters"))
})

test_that("kmer_analysis threshold parameter works correctly", {
  skip_on_cran()
  
  motif_matrix <- matrix(
    c(10, 0, 0,
      2, 0, 0,    # Below threshold of 5
      0, 8, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      c("motif_A", "motif_B", "motif_C"),
      c("Cluster_1", "Cluster_2", "Cluster_3")
    )
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:3,
    dominant_class = c("Class_A", "Class_B", "Class_A"),
    stringsAsFactors = FALSE
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  # With threshold = 5
  result <- kmer_analysis(
    cluster_result = cluster_result,
    motif_matrix = motif_matrix,
    threshold = 5
  )
  
  # motif_B should not be considered present
  unique_motifs <- result$unique_cluster_motifs$motif
  expect_true("motif_A" %in% unique_motifs)
  expect_false("motif_B" %in% unique_motifs)
  expect_true("motif_C" %in% unique_motifs)
})

test_that("kmer_analysis handles validation data correctly", {
  skip_on_cran()
  
  motif_matrix <- matrix(
    c(10, 0, 5, 0, 8, 0),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      c("motif_A", "motif_B"),
      c("Cluster_1", "Cluster_2", "Cluster_3")
    )
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:3,
    dominant_class = c("Class_A", "Class_B", "Class_A"),
    stringsAsFactors = FALSE
  )
  
  # Add validation data with proper structure
  # Includes data_result with kmers for training rank calculation
  training_data <- data.frame(
    motif_A = c(10, 5, 3),
    motif_B = c(8, 0, 7),
    CLASS = c("Class_A", "Class_B", "Class_A"),
    stringsAsFactors = FALSE
  )
  
  validation_data <- data.frame(
    motif_A = c(5, 0, 3),
    motif_B = c(0, 7, 0),
    CLASS = c("Class_A", "Class_B", "Class_A"),
    stringsAsFactors = FALSE
  )
  
  cluster_result <- list(
    cluster_summary = cluster_summary,
    data_result = list(kmers = training_data),  # Added training data
    classification_result = list(validation_data = validation_data)
  )
  
  result <- kmer_analysis(
    cluster_result = cluster_result,
    motif_matrix = motif_matrix,
    threshold = 0
  )
  
  # Check validation flag
  expect_true(result$has_external_validation)
  expect_type(result$validation, "list")
  
  # Check validation components
  expect_named(result$validation, c(
    "unique_class_motifs",
    "class_frequency_ranking",
    "class_frequency_matrix",
    "motifs_by_class_rank",
    "n_sequences",
    "n_classes",
    "class_names"
  ))
  
  # Check validation statistics
  expect_equal(result$validation$n_sequences, 3)
  expect_equal(result$validation$n_classes, 2)
  expect_setequal(result$validation$class_names, c("Class_A", "Class_B"))
  
  # Check that motifs_by_class_rank was created for validation
  expect_s3_class(result$validation$motifs_by_class_rank, "data.frame")
  expect_true("motif" %in% colnames(result$validation$motifs_by_class_rank))
})

test_that("kmer_analysis handles empty results gracefully", {
  skip_on_cran()
  
  # Matrix where all motifs appear in all clusters
  motif_matrix <- matrix(
    c(5, 5, 5,
      3, 3, 3),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      c("motif_A", "motif_B"),
      c("Cluster_1", "Cluster_2", "Cluster_3")
    )
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:3,
    dominant_class = rep("Class_A", 3),  # All same class
    stringsAsFactors = FALSE
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  result <- kmer_analysis(
    cluster_result = cluster_result,
    motif_matrix = motif_matrix,
    threshold = 0
  )
  
  # Should have no unique CLUSTER motifs (all appear in multiple clusters)
  expect_equal(nrow(result$unique_cluster_motifs), 0)
  
  # Should have unique CLASS motifs (all motifs only in Class_A)
  expect_equal(nrow(result$unique_class_motifs), 2)
  expect_true(all(result$unique_class_motifs$class == "Class_A"))
  
  # But should have frequency rankings
  expect_equal(nrow(result$cluster_frequency_ranking), 2)
  expect_true("n_clusters" %in% colnames(result$cluster_frequency_ranking))
})

test_that("kmer_analysis validates input structure", {
  skip_on_cran()
  
  motif_matrix <- matrix(1:9, nrow = 3)
  
  # Missing cluster_summary
  cluster_result <- list()
  
  expect_error(
    kmer_analysis(
      cluster_result = cluster_result,
      motif_matrix = motif_matrix
    )
  )
  
  # Invalid motif_matrix (not a matrix)
  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = 1:3,
      dominant_class = c("A", "B", "C")
    )
  )
  
  expect_error(
    kmer_analysis(
      cluster_result = cluster_result,
      motif_matrix = "not_a_matrix"
    )
  )
})

test_that("kmer_analysis cluster_to_class mapping is correct", {
  skip_on_cran()
  
  motif_matrix <- matrix(
    1:6, nrow = 2,
    dimnames = list(
      c("motif_A", "motif_B"),
      c("Cluster_1", "Cluster_2", "Cluster_3")
    )
  )
  
  cluster_summary <- data.frame(
    cluster_id = 1:3,
    dominant_class = c("Class_X", "Class_Y", "Class_X"),
    stringsAsFactors = FALSE
  )
  
  cluster_result <- list(cluster_summary = cluster_summary)
  
  result <- kmer_analysis(
    cluster_result = cluster_result,
    motif_matrix = motif_matrix,
    threshold = 0
  )
  
  # Check mapping structure
  expect_type(result$cluster_to_class, "character")
  expect_named(result$cluster_to_class, 
               c("Cluster_1", "Cluster_2", "Cluster_3"))
  expect_equal(result$cluster_to_class[["Cluster_1"]], "Class_X")
  expect_equal(result$cluster_to_class[["Cluster_2"]], "Class_Y")
  expect_equal(result$cluster_to_class[["Cluster_3"]], "Class_X")
})
