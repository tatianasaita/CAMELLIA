# tests/testthat/test-select_sample_train_validation.R

library(testthat)

test_that("select_sample_train_validation validates inputs", {
  # Setup minimal cluster_result
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = paste0("seq", 1:20),
        length = rep(c(500, 1000), each = 10),
        class = rep(c("ClassA", "ClassB"), each = 10),
        stringsAsFactors = FALSE
      ),
      kmers = matrix(runif(20 * 10), nrow = 20, ncol = 10,
                     dimnames = list(paste0("seq", 1:20), paste0("kmer", 1:10)))
    )
  )
  
  # Test: min_size filter removes all sequences
  expect_error(
    select_sample_train_validation(cluster_result, k_per_class = 5, min_size = 2000, verbose = FALSE),
    "No sequences available with length >= 2000"
  )
  
  # Test: k_per_class larger than available sequences (should adjust automatically)
  expect_no_error({
    result <- select_sample_train_validation(cluster_result, k_per_class = 100, min_size = 500, verbose = FALSE)
  })
  
  # Optional: verify it selected all available sequences
  result <- select_sample_train_validation(cluster_result, k_per_class = 100, min_size = 500, verbose = FALSE)
  expect_true(nrow(result$classification_metadata) <= 20)
})

test_that("select_sample_train_validation returns correct structure", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = paste0("seq", 1:20),
        length = rep(1000, 20),
        class = rep(c("ClassA", "ClassB"), each = 10),
        stringsAsFactors = FALSE
      ),
      kmers = matrix(runif(20 * 10), nrow = 20, ncol = 10,
                     dimnames = list(paste0("seq", 1:20), paste0("kmer", 1:10)))
    )
  )
  
  result <- select_sample_train_validation(cluster_result, k_per_class = 3, min_size = 800, verbose = FALSE)
  
  # Test class
  expect_s3_class(result, "sample_selection")
  expect_s3_class(result, "list")
  
  # Test structure
  expect_named(result, c("classification_dataset", "validation_dataset", 
                         "classification_metadata", "validation_metadata",
                         "true_labels_classification", "true_labels_validation"))
  
  # Test attributes
  expect_true(is.numeric(attr(result, "n_classification")))
  expect_true(is.numeric(attr(result, "n_validation")))
  expect_equal(attr(result, "validation_type"), "internal")
  expect_equal(attr(result, "min_size"), 800)
  expect_equal(attr(result, "k_per_class"), 3)
})

test_that("select_sample_train_validation splits data correctly", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = paste0("seq", 1:20),
        length = rep(1000, 20),
        class = rep(c("ClassA", "ClassB"), each = 10),
        stringsAsFactors = FALSE
      ),
      kmers = matrix(runif(20 * 10), nrow = 20, ncol = 10,
                     dimnames = list(paste0("seq", 1:20), paste0("kmer", 1:10)))
    )
  )
  
  result <- select_sample_train_validation(cluster_result, k_per_class = 3, verbose = FALSE)
  
  # Test no overlap between classification and validation
  overlap <- intersect(
    result$classification_metadata$sequence_name,
    result$validation_metadata$sequence_name
  )
  expect_length(overlap, 0)
  
  # Test classification has correct number per class
  expect_equal(nrow(result$classification_metadata), 6)  # 3 per class * 2 classes
  
  # Test validation has remaining sequences
  expect_equal(nrow(result$validation_metadata), 14)  # 20 - 6
  
  # Test k-mer dimensions match metadata
  expect_equal(nrow(result$classification_dataset), nrow(result$classification_metadata))
  expect_equal(nrow(result$validation_dataset), nrow(result$validation_metadata))
})

test_that("select_sample_train_validation handles min_size filter", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = paste0("seq", 1:20),
        length = c(rep(500, 10), rep(1000, 10)),
        class = rep(c("ClassA", "ClassB"), each = 10),
        stringsAsFactors = FALSE
      ),
      kmers = matrix(runif(20 * 10), nrow = 20, ncol = 10,
                     dimnames = list(paste0("seq", 1:20), paste0("kmer", 1:10)))
    )
  )
  
  result <- select_sample_train_validation(cluster_result, k_per_class = 3, min_size = 800, verbose = FALSE)
  
  # Test: only sequences >= min_size in classification
  expect_true(all(result$classification_metadata$length >= 800))
  
  # Test: validation includes sequences < min_size
  expect_true(any(result$validation_metadata$length < 800))
})

test_that("select_sample_train_validation verbose parameter works", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = paste0("seq", 1:10),
        length = rep(1000, 10),
        class = rep("ClassA", 10),
        stringsAsFactors = FALSE
      ),
      kmers = matrix(runif(10 * 5), nrow = 10, ncol = 5,
                     dimnames = list(paste0("seq", 1:10), paste0("kmer", 1:5)))
    )
  )
  
  # Test: verbose = FALSE produces no output
  expect_silent({
    result <- select_sample_train_validation(cluster_result, k_per_class = 3, verbose = FALSE)
  })
  
  # Test: verbose = TRUE produces output
  expect_output(
    select_sample_train_validation(cluster_result, k_per_class = 3, verbose = TRUE),
    "Sequences with length"
  )
})

test_that("select_sample_train_validation handles empty validation set", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = paste0("seq", 1:6),
        length = rep(1000, 6),
        class = rep(c("ClassA", "ClassB"), each = 3),
        stringsAsFactors = FALSE
      ),
      kmers = matrix(runif(6 * 5), nrow = 6, ncol = 5,
                     dimnames = list(paste0("seq", 1:6), paste0("kmer", 1:5)))
    )
  )
  
  # Select all sequences for classification
  expect_warning(
    result <- select_sample_train_validation(cluster_result, k_per_class = 3, verbose = FALSE),
    "No sequences available for internal validation"
  )
  
  expect_equal(nrow(result$validation_metadata), 0)
  expect_equal(nrow(result$validation_dataset), 0)
})

test_that("select_sample_train_validation saves files", {
  # Use temporary directory for tests
  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = paste0("seq", 1:10),
        length = rep(1000, 10),
        class = rep("ClassA", 10),
        stringsAsFactors = FALSE
      ),
      kmers = matrix(runif(10 * 5), nrow = 10, ncol = 5,
                     dimnames = list(paste0("seq", 1:10), paste0("kmer", 1:5)))
    )
  )
  
  result <- select_sample_train_validation(cluster_result, k_per_class = 3, verbose = FALSE)
  
  # Test: files are created
  expect_true(file.exists("true_labels_classification.RData"))
  expect_true(file.exists("classification_dataset.RData"))
  expect_true(file.exists("true_labels_validation.RData"))
  expect_true(file.exists("validation_dataset.RData"))
  expect_true(file.exists("backup_checkpoint1.RData"))
  
  # Clean up
  file.remove("true_labels_classification.RData")
  file.remove("classification_dataset.RData")
  file.remove("true_labels_validation.RData")
  file.remove("validation_dataset.RData")
  file.remove("backup_checkpoint1.RData")
})

