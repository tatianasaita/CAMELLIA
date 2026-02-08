# tests/testthat/test-create_dendrogram.R

library("testthat")
test_that("create_dendrogram returns correct structure", {
  # Create test data
  data <- data.frame(
    kmer1 = rnorm(50),
    kmer2 = rnorm(50),
    kmer3 = rnorm(50),
    CLASS = rep(c("A", "B"), each = 25)
  )
  
  result <- create_dendrogram(data)
  
  # Check class
  expect_s3_class(result, "dendrogram_result")
  expect_type(result, "list")
  
  # Check required components
  expect_named(result, c("dendrogram", "hclust", "order", "labels", 
                         "sequence_names", "colors", "base_colors"))
  
  # Check component types
  expect_s3_class(result$dendrogram, "dendrogram")
  expect_s3_class(result$hclust, "hclust")
  expect_type(result$order, "integer")
  expect_type(result$labels, "character")
  expect_type(result$colors, "character")
  expect_type(result$base_colors, "character")
  
  # Check dimensions
  expect_length(result$order, nrow(data))
  expect_length(result$labels, nrow(data))
  expect_length(result$colors, nrow(data))
  expect_length(result$base_colors, 2)  # Two classes
})

test_that("create_dendrogram handles sequence_names correctly", {
  data <- data.frame(
    kmer1 = rnorm(30),
    kmer2 = rnorm(30),
    CLASS = rep(c("A", "B", "C"), each = 10)
  )
  
  seq_names <- paste0("seq_", 1:30)
  result <- create_dendrogram(data, sequence_names = seq_names)
  
  expect_type(result$sequence_names, "character")
  expect_length(result$sequence_names, 30)
  expect_true(all(result$sequence_names %in% seq_names))
})

test_that("create_dendrogram works with different distance methods", {
  data <- data.frame(
    kmer1 = rnorm(20),
    kmer2 = rnorm(20),
    CLASS = rep(c("A", "B"), each = 10)
  )
  
  methods <- c("euclidean", "manhattan", "maximum")
  
  for (method in methods) {
    result <- create_dendrogram(data, dist_method = method)
    expect_s3_class(result, "dendrogram_result")
  }
})

test_that("create_dendrogram works with different clustering methods", {
  data <- data.frame(
    kmer1 = rnorm(20),
    kmer2 = rnorm(20),
    CLASS = rep(c("A", "B"), each = 10)
  )
  
  methods <- c("ward.D2", "complete", "average", "single")
  
  for (method in methods) {
    result <- create_dendrogram(data, hclust_method = method)
    expect_s3_class(result, "dendrogram_result")
  }
})

test_that("create_dendrogram handles multiple classes correctly", {
  # Test with 2 classes
  data_2 <- data.frame(
    kmer1 = rnorm(20),
    kmer2 = rnorm(20),
    CLASS = rep(c("A", "B"), each = 10)
  )
  result_2 <- create_dendrogram(data_2)
  expect_length(result_2$base_colors, 2)
  
  # Test with 5 classes (uses palette.colors)
  data_5 <- data.frame(
    kmer1 = rnorm(50),
    kmer2 = rnorm(50),
    CLASS = rep(c("A", "B", "C", "D", "E"), each = 10)
  )
  result_5 <- create_dendrogram(data_5)
  expect_length(result_5$base_colors, 5)
  
  # Test with 15 classes (>12, should use rainbow)
  data_15 <- data.frame(
    kmer1 = rnorm(150),
    kmer2 = rnorm(150),
    CLASS = rep(LETTERS[1:15], each = 10)
  )
  result_15 <- create_dendrogram(data_15)
  expect_length(result_15$base_colors, 15)
})

test_that("create_dendrogram saves output file correctly", {
  skip_on_cran()
  
  data <- data.frame(
    kmer1 = rnorm(20),
    kmer2 = rnorm(20),
    CLASS = rep(c("A", "B"), each = 10)
  )
  
  temp_file <- tempfile(fileext = ".png")
  
  result <- create_dendrogram(data, output = temp_file)
  
  expect_true(file.exists(temp_file))
  expect_gt(file.size(temp_file), 0)
  
  # Clean up
  unlink(temp_file)
})

test_that("create_dendrogram color assignment is consistent", {
  data <- data.frame(
    kmer1 = rnorm(30),
    kmer2 = rnorm(30),
    CLASS = rep(c("A", "B", "C"), each = 10)
  )
  
  result <- create_dendrogram(data)
  
  # Check that base_colors are named
  expect_true(!is.null(names(result$base_colors)))
  expect_setequal(names(result$base_colors), c("A", "B", "C"))
  
  # Check that colors match classes
  expect_length(result$colors, 30)
  
  # Check color consistency
  for (i in seq_along(result$labels)) {
    expected_color <- result$base_colors[result$labels[i]]
    expect_equal(result$colors[i], expected_color, ignore_attr = TRUE)
  }
})

test_that("print method works correctly", {
  data <- data.frame(
    kmer1 = rnorm(20),
    kmer2 = rnorm(20),
    CLASS = rep(c("A", "B"), each = 10)
  )
  
  result <- create_dendrogram(data)
  
  # Capture output
  output <- capture.output(print(result))
  
  # Check key elements are present
  expect_true(any(grepl("Dendrogram Analysis Result", output)))
  expect_true(any(grepl("Number of samples:", output)))
  expect_true(any(grepl("Number of classes:", output)))
  expect_true(any(grepl("Tree height:", output)))
  expect_true(any(grepl("Color mapping by class:", output)))
})

test_that("create_dendrogram order is valid", {
  data <- data.frame(
    kmer1 = rnorm(25),
    kmer2 = rnorm(25),
    CLASS = rep(c("A", "B"), c(10, 15))
  )
  
  result <- create_dendrogram(data)
  
  # Check order contains all indices
  expect_setequal(result$order, 1:25)
  
  # Check no duplicates
  expect_equal(length(result$order), length(unique(result$order)))
  
  # Check order is valid permutation
  expect_true(all(result$order >= 1 & result$order <= 25))
})

test_that("create_dendrogram handles edge cases", {
  # Minimum valid data (2 rows)
  data_min <- data.frame(
    kmer1 = c(1, 2),
    kmer2 = c(3, 4),
    CLASS = c("A", "B")
  )
  
  result <- create_dendrogram(data_min)
  expect_s3_class(result, "dendrogram_result")
  expect_length(result$order, 2)
})

test_that("create_dendrogram works without sequence_names", {
  data <- data.frame(
    kmer1 = rnorm(15),
    kmer2 = rnorm(15),
    CLASS = rep(c("A", "B", "C"), each = 5)
  )
  
  result <- create_dendrogram(data)
  
  # Should create default sequence names
  expect_type(result$sequence_names, "character")
  expect_length(result$sequence_names, 15)
})

