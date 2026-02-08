# tests/testthat/test-cluster_dendrogram.R

library(testthat)

setup_test_dendrogram <- function() {
  dend <- as.dendrogram(hclust(dist(1:10)))
  return(dend)
}


test_that("cluster_dendrogram returns correct structure", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE
  )
  
  expect_s3_class(result, "cluster_dendrogram_result")
  expect_type(result, "list")
  expect_named(result, c(
    "dendrogram", "clusters", "cluster_summary", "element_assignment",
    "data_result", "cluster_assignment_dendro_order", "unassigned_elements",
    "n_unassigned", "hom_thresh", "min_size", "n_elements", "valid_clusters"
  ))
})


test_that("cluster_dendrogram handles homogeneity threshold correctly", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE
  )
  
  expect_true(all(result$cluster_summary$homogeneity >= 0.7))
  expect_equal(result$hom_thresh, 0.7)
})


test_that("cluster_dendrogram respects minimum cluster size", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  min_size <- 3
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = min_size,
    verbose = FALSE
  )
  
  expect_true(all(result$cluster_summary$n_elements >= min_size))
  expect_equal(result$min_size, min_size)
})


test_that("cluster_dendrogram handles rare classes correctly", {
  dend <- as.dendrogram(hclust(dist(1:8)))
  class_labels <- c("A", "A", "A", "B", "B", "B", "D", "D")
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE
  )
  
  d_clusters <- result$cluster_summary[result$cluster_summary$dominant_class == "D", ]
  expect_true(nrow(d_clusters) > 0)
})


test_that("cluster_dendrogram assigns all elements", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE
  )
  
  expect_equal(result$n_unassigned, 0)
  expect_equal(length(result$unassigned_elements), 0)
  expect_true(all(!is.na(result$cluster_assignment_dendro_order)))
})


test_that("cluster_dendrogram creates valid cluster_summary", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE
  )
  
  expect_s3_class(result$cluster_summary, "data.frame")
  expect_named(result$cluster_summary, c(
    "cluster_id", "n_elements", "dominant_class", "homogeneity",
    "n_classes", "class_composition", "is_complete_class"
  ))
  expect_true(all(result$cluster_summary$homogeneity >= 0 & 
                    result$cluster_summary$homogeneity <= 1))
  expect_true(all(result$cluster_summary$n_elements > 0))
})


test_that("cluster_dendrogram creates valid element_assignment", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  seq_names <- paste0("seq", 1:10)
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE,
    sequence_names = seq_names
  )
  
  expect_s3_class(result$element_assignment, "data.frame")
  expect_equal(nrow(result$element_assignment), length(class_labels))
  expect_named(result$element_assignment, c(
    "sequence_name", "dendro_index", "class", "cluster"
  ))
  expect_equal(result$element_assignment$sequence_name, seq_names)
})


test_that("cluster_dendrogram handles single class correctly", {
  dend <- setup_test_dendrogram()
  class_labels <- rep("A", 10)
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE
  )
  
  expect_equal(nrow(result$cluster_summary), 1)
  expect_equal(result$cluster_summary$homogeneity, 1.0)
  expect_equal(result$cluster_summary$dominant_class, "A")
})


test_that("cluster_dendrogram handles data_result parameter", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  
  data_result <- list(
    kmers = matrix(rnorm(100), ncol = 10),
    metadata = data.frame(id = 1:10, class = class_labels)
  )
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE,
    data_result = data_result
  )
  
  expect_true("cluster" %in% names(result$data_result$metadata))
  expect_equal(length(result$data_result$metadata$cluster), 10)
})


test_that("cluster_dendrogram handles verbose parameter", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  
  expect_silent(
    cluster_dendrogram(
      dendrogram = dend,
      class_labels = class_labels,
      hom_thresh = 0.7,
      min_size = 3,
      verbose = FALSE
    )
  )
  
  expect_error(
    cluster_dendrogram(
      dendrogram = dend,
      class_labels = class_labels,
      hom_thresh = 0.7,
      min_size = 3,
      verbose = TRUE
    ),
    NA
  )
})


test_that("cluster_dendrogram handles minimum valid input", {
  dend <- as.dendrogram(hclust(dist(1:3)))
  class_labels <- c("A", "A", "A")
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE
  )
  
  expect_s3_class(result, "cluster_dendrogram_result")
  expect_equal(result$n_elements, 3)
})


test_that("cluster_dendrogram preserves dendrogram order", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE
  )
  
  all_indices <- unlist(lapply(result$clusters, function(cl) cl$indices))
  expect_true(all(all_indices >= 1 & all_indices <= length(class_labels)))
  expect_equal(length(unique(all_indices)), length(all_indices))
})


test_that("cluster_dendrogram handles NULL sequence_names", {
  dend <- setup_test_dendrogram()
  class_labels <- c("A", "A", "A", "B", "B", "B", "C", "C", "C", "C")
  
  result <- cluster_dendrogram(
    dendrogram = dend,
    class_labels = class_labels,
    hom_thresh = 0.7,
    min_size = 3,
    verbose = FALSE,
    sequence_names = NULL
  )
  
  expect_s3_class(result$element_assignment, "data.frame")
  expect_true(!is.null(result$element_assignment$sequence_name))
})
