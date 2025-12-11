# tests/testthat/test-calculate_cluster_motifs.R
#
# Essential tests for calculate_cluster_motifs function following CRAN standards
# These tests cover:
# - Input validation (object class, data structure)
# - Basic functionality (aggregation and normalization)
# - Output structure and format
# - S3 methods (print, summary)
# - Edge cases (zero-range columns)

test_that("calculate_cluster_motifs validates inputs correctly", {
  # Test invalid object class
  expect_error(
    calculate_cluster_motifs(list()),
    "'x' must be an object of class 'cluster_dendrogram_result'"
  )

  # Test missing data_result
  x <- structure(list(), class = "cluster_dendrogram_result")
  expect_error(
    calculate_cluster_motifs(x),
    "'x' must contain 'data_result' element"
  )

  # Test missing cluster column
  x <- structure(
    list(
      data_result = list(
        kmers = data.frame(AAA = 1:5, TTT = 2:6),
        metadata = data.frame(id = 1:5)
      )
    ),
    class = "cluster_dendrogram_result"
  )
  expect_error(
    calculate_cluster_motifs(x),
    "'data_result\\$metadata' must contain a 'cluster' column"
  )
})

test_that("calculate_cluster_motifs aggregates k-mers correctly", {
  # Create simple test data
  kmers_df <- data.frame(
    AAA = c(10, 20, 30, 40),
    TTT = c(5, 10, 15, 20),
    CLASS = c("A", "A", "B", "B")
  )

  metadata <- data.frame(
    id = 1:4,
    cluster = c(1, 1, 2, 2)
  )

  x <- structure(
    list(
      data_result = list(
        kmers = kmers_df,
        metadata = metadata
      )
    ),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(x)

  # Check basic structure
  expect_s3_class(result, "cluster_motifs")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)  # 2 k-mers (AAA, TTT)
  expect_equal(ncol(result), 2)  # 2 clusters

  # Check column names are ordered
  expect_equal(colnames(result), c("Cluster_1", "Cluster_2"))

  # Check row names
  expect_equal(rownames(result), c("AAA", "TTT"))
})

test_that("calculate_cluster_motifs normalizes correctly", {
  # Create data with known values
  kmers_df <- data.frame(
    AAA = c(0, 100, 0, 100),
    TTT = c(50, 50, 50, 50),
    CLASS = c("A", "A", "B", "B")
  )

  metadata <- data.frame(
    id = 1:4,
    cluster = c(1, 1, 2, 2)
  )

  x <- structure(
    list(
      data_result = list(
        kmers = kmers_df,
        metadata = metadata
      )
    ),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(x)

  # All values should be between 0 and 1
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("calculate_cluster_motifs handles zero-range columns", {
  # Create data where all k-mer sums in a cluster are identical
  kmers_df <- data.frame(
    AAA = c(10, 10, 10, 10),
    TTT = c(10, 10, 10, 10),
    CLASS = c("A", "A", "B", "B")
  )

  metadata <- data.frame(
    id = 1:4,
    cluster = c(1, 1, 2, 2)
  )

  x <- structure(
    list(
      data_result = list(
        kmers = kmers_df,
        metadata = metadata
      )
    ),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(x)

  # Both clusters have identical sums for all k-mers
  # Cluster 1: AAA=20, TTT=20 (range=0)
  # Cluster 2: AAA=20, TTT=20 (range=0)
  # Should all be set to 0.5
  expect_equal(result["AAA", "Cluster_1"], 0.5)
  expect_equal(result["TTT", "Cluster_1"], 0.5)
  expect_equal(result["AAA", "Cluster_2"], 0.5)
  expect_equal(result["TTT", "Cluster_2"], 0.5)

  # No NA or NaN values
  expect_false(any(is.na(result)))
  expect_false(any(is.nan(as.matrix(result))))
})

test_that("calculate_cluster_motifs normalizes with different ranges", {
  # Create data with different k-mer abundances
  kmers_df <- data.frame(
    AAA = c(10, 10, 20, 20),
    TTT = c(5, 5, 5, 5),
    CLASS = c("A", "A", "B", "B")
  )

  metadata <- data.frame(
    id = 1:4,
    cluster = c(1, 1, 2, 2)
  )

  x <- structure(
    list(
      data_result = list(
        kmers = kmers_df,
        metadata = metadata
      )
    ),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(x)

  # Cluster 1: AAA=20, TTT=10 → normalized: AAA=1.0, TTT=0.0
  expect_equal(result["AAA", "Cluster_1"], 1.0)
  expect_equal(result["TTT", "Cluster_1"], 0.0)

  # Cluster 2: AAA=40, TTT=10 → normalized: AAA=1.0, TTT=0.0
  expect_equal(result["AAA", "Cluster_2"], 1.0)
  expect_equal(result["TTT", "Cluster_2"], 0.0)

  # All values should be between 0 and 1
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("calculate_cluster_motifs removes CLASS column", {
  # Data with CLASS column
  kmers_df <- data.frame(
    CLASS = c("A", "A", "B", "B"),
    AAA = c(10, 20, 30, 40),
    TTT = c(5, 10, 15, 20)
  )

  metadata <- data.frame(
    id = 1:4,
    cluster = c(1, 1, 2, 2)
  )

  x <- structure(
    list(
      data_result = list(
        kmers = kmers_df,
        metadata = metadata
      )
    ),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(x)

  # Should only have k-mer columns
  expect_equal(rownames(result), c("AAA", "TTT"))
  expect_false("CLASS" %in% rownames(result))
})

test_that("calculate_cluster_motifs accepts external data_result", {
  # Create objects
  x <- structure(list(), class = "cluster_dendrogram_result")

  data_result <- list(
    kmers = data.frame(
      AAA = c(10, 20, 30),
      TTT = c(5, 10, 15)
    ),
    metadata = data.frame(
      id = 1:3,
      cluster = c(1, 2, 2)
    )
  )

  result <- calculate_cluster_motifs(x, data_result = data_result)

  expect_s3_class(result, "cluster_motifs")
  expect_equal(ncol(result), 2)  # 2 clusters
})

test_that("print method works correctly", {
  kmers_df <- data.frame(
    AAA = c(10, 20),
    TTT = c(5, 10)
  )

  metadata <- data.frame(
    id = 1:2,
    cluster = c(1, 1)
  )

  x <- structure(
    list(
      data_result = list(
        kmers = kmers_df,
        metadata = metadata
      )
    ),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(x)

  expect_output(print(result), "Cluster Motif Matrix")
  expect_output(print(result), "Dimensions:")
  expect_output(print(result), "Normalization: Min-Max")
})

test_that("summary method works correctly", {
  kmers_df <- data.frame(
    AAA = c(10, 20, 30),
    TTT = c(5, 10, 15)
  )

  metadata <- data.frame(
    id = 1:3,
    cluster = c(1, 2, 2)
  )

  x <- structure(
    list(
      data_result = list(
        kmers = kmers_df,
        metadata = metadata
      )
    ),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(x)

  expect_output(summary(result), "Cluster Motif Summary")
  expect_output(summary(result), "Matrix dimensions:")
  expect_output(summary(result), "Per-cluster statistics:")
  expect_output(summary(result), "min=")
  expect_output(summary(result), "median=")
  expect_output(summary(result), "mean=")
  expect_output(summary(result), "max=")
})

test_that("calculate_cluster_motifs handles single cluster", {
  kmers_df <- data.frame(
    AAA = c(10, 20, 30),
    TTT = c(5, 10, 15)
  )

  metadata <- data.frame(
    id = 1:3,
    cluster = c(1, 1, 1)
  )

  x <- structure(
    list(
      data_result = list(
        kmers = kmers_df,
        metadata = metadata
      )
    ),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(x)

  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "Cluster_1")
})

test_that("calculate_cluster_motifs preserves cluster order", {
  kmers_df <- data.frame(
    AAA = c(10, 20, 30, 40),
    TTT = c(5, 10, 15, 20)
  )

  # Non-sequential cluster IDs
  metadata <- data.frame(
    id = 1:4,
    cluster = c(3, 1, 3, 1)
  )

  x <- structure(
    list(
      data_result = list(
        kmers = kmers_df,
        metadata = metadata
      )
    ),
    class = "cluster_dendrogram_result"
  )

  result <- calculate_cluster_motifs(x)

  # Should be sorted: Cluster_1, Cluster_3
  expect_equal(colnames(result), c("Cluster_1", "Cluster_3"))
})

