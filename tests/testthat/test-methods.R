# tests/testthat/test-methods.R

library(testthat)

# Helper function to capture output
capture_print <- function(x) {
  paste(capture.output(print(x)), collapse = "\n")
}


# Tests for print.kmer_data

test_that("print.kmer_data displays correct structure", {
  skip_if_not_installed("Biostrings")

  # Create mock kmer_data object
  kmer_data <- list(
    kmers = data.frame(
      AT = c(5, 3, 2),
      CG = c(1, 4, 6),
      CLASS = c("A", "A", "B"),
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      sequence_name = c("seq1", "seq2", "seq3"),
      length = c(100, 150, 120),
      class = c("A", "A", "B"),
      stringsAsFactors = FALSE
    )
  )
  class(kmer_data) <- c("kmer_data", "list")

  output <- capture_print(kmer_data)

  expect_match(output, "K-mer Data Object")
  expect_match(output, "K-mer Matrix:")
  expect_match(output, "3 sequences x 2 k-mers")
  expect_match(output, "Metadata:")
  expect_match(output, "Class Distribution:")
  expect_match(output, "A: 2 sequences")
  expect_match(output, "B: 1 sequences")
})

test_that("print.kmer_data returns invisibly", {
  kmer_data <- list(
    kmers = data.frame(AT = 1, CLASS = "A"),
    metadata = data.frame(sequence_name = "seq1", length = 100, class = "A")
  )
  class(kmer_data) <- c("kmer_data", "list")

  result <- expect_invisible(print(kmer_data))
  expect_s3_class(result, "kmer_data")
})

# Tests for print.dendrogram_result

test_that("print.dendrogram_result displays correct information", {
  # Create mock dendrogram_result
  hc <- hclust(dist(matrix(1:6, ncol = 2)))
  dend <- as.dendrogram(hc)

  dend_result <- list(
    dendrogram = dend,
    hclust = hc,
    order = c(1, 2, 3),
    labels = c("A", "B", "A"),
    sequence_names = c("seq1", "seq2", "seq3"),
    colors = c("red", "blue", "red"),
    base_colors = c(A = "red", B = "blue")
  )
  class(dend_result) <- c("dendrogram_result", "list")

  output <- capture_print(dend_result)

  expect_match(output, "Dendrogram Analysis Result")
  expect_match(output, "Number of samples: 3")
  expect_match(output, "Number of classes: 2")
  expect_match(output, "Tree height:")
  expect_match(output, "Color mapping by class:")
})

test_that("print.dendrogram_result handles missing sequence_names", {
  hc <- hclust(dist(matrix(1:6, ncol = 2)))
  dend <- as.dendrogram(hc)

  dend_result <- list(
    dendrogram = dend,
    hclust = hc,
    order = c(1, 2, 3),
    labels = c("A", "B", "A"),
    sequence_names = NULL,
    colors = c("red", "blue", "red"),
    base_colors = c(A = "red", B = "blue")
  )
  class(dend_result) <- c("dendrogram_result", "list")

  output <- capture_print(dend_result)

  expect_false(grepl("sequence_names:", output, fixed = TRUE))
})

# Tests for print.cluster_dendrogram_result

test_that("print.cluster_dendrogram_result displays summary correctly", {
  cluster_result <- list(
    n_elements = 10,
    cluster_summary = data.frame(
      cluster_id = c(1, 2),
      size = c(5, 3),
      homogeneity = c(1.0, 0.8),
      dominant_class = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    n_unassigned = 2,
    min_size = 3,
    hom_thresh = 0.7
  )
  class(cluster_result) <- c("cluster_dendrogram_result", "list")

  output <- capture_print(cluster_result)

  expect_match(output, "Dendrogram Clustering Result")
  expect_match(output, "Total elements:\\s+10")
  expect_match(output, "Total clusters:\\s+2")
  expect_match(output, "Unassigned:\\s+2")
  expect_match(output, "Min size:\\s+3")
  expect_match(output, "Hom threshold:\\s+0.700")
  expect_match(output, "Cluster Summary:")
  expect_match(output, "WARNING.*could not be assigned")
})

test_that("print.cluster_dendrogram_result handles no unassigned", {
  cluster_result <- list(
    n_elements = 10,
    cluster_summary = data.frame(cluster_id = 1, size = 10,
                                 homogeneity = 1.0, dominant_class = "A"),
    n_unassigned = 0,
    min_size = 3,
    hom_thresh = 0.7
  )
  class(cluster_result) <- c("cluster_dendrogram_result", "list")

  output <- capture_print(cluster_result)

  expect_false(grepl("WARNING", output))
})

# Tests for print.cluster_motifs_result

test_that("print.cluster_motifs_result displays matrix info", {
  motif_matrix <- data.frame(
    cluster1 = c(0.5, 0.8, 0.2),
    cluster2 = c(0.9, 0.1, 0.6),
    row.names = c("ATCG", "GCTA", "TACG")
  )
  class(motif_matrix) <- c("cluster_motifs_result", "data.frame")

  output <- capture_print(motif_matrix)

  expect_match(output, "Cluster Motif Matrix")
  expect_match(output, "Dimensions: 3 k-mers x 2 clusters")
  expect_match(output, "Normalization: Min-Max")
  expect_match(output, "Value range: \\[0, 1\\]")
})

# Tests for print.select_motifs

test_that("print.select_motifs displays selection summary", {
  selected <- list(
    A = c("ATCG", "GCTA"),
    B = c("TACG", "CGAT", "ATAT")
  )
  class(selected) <- "select_motifs"
  attr(selected, "n") <- 5
  attr(selected, "k") <- 2
  attr(selected, "m") <- 4
  attr(selected, "case") <- "n < m: select top motifs per cluster"

  output <- capture_print(selected)

  expect_match(output, "Selected Motifs Summary")
  expect_match(output, "Total motifs \\(n\\):\\s+5")
  expect_match(output, "Classes \\(k\\):\\s+2")
  expect_match(output, "Clusters \\(m\\):\\s+4")
  expect_match(output, "Motifs per class:")
  expect_match(output, "A: 2")
  expect_match(output, "B: 3")
})

test_that("print.select_motifs handles empty results", {
  selected <- list()
  class(selected) <- c("select_motifs", "empty_result")
  attr(selected, "reason") <- "No motifs available"
  attr(selected, "n") <- 0
  attr(selected, "k") <- 0
  attr(selected, "m") <- 0
  attr(selected, "case") <- "empty"

  output <- capture_print(selected)

  expect_match(output, "Result: Empty")
  expect_match(output, "No motifs available")
})

# Tests for print.sample_selection

test_that("print.sample_selection displays sample info", {
  sample_sel <- list(
    classification_metadata = data.frame(
      sequence_name = paste0("seq", 1:10),
      class = rep(c("A", "B"), each = 5)
    ),
    validation_metadata = data.frame(
      sequence_name = paste0("val", 1:4),
      class = rep(c("A", "B"), each = 2)
    )
  )
  class(sample_sel) <- "sample_selection"
  attr(sample_sel, "n_classification") <- 10
  attr(sample_sel, "n_validation") <- 4
  attr(sample_sel, "validation_type") <- "internal"
  attr(sample_sel, "min_size") <- 50
  attr(sample_sel, "seq_per_class") <- 5

  output <- capture_print(sample_sel)

  expect_match(output, "Sample Selection Summary")
  expect_match(output, "Classification sequences: 10")
  expect_match(output, "Validation sequences:\\s+4")
  expect_match(output, "Validation type:\\s+internal")
  expect_match(output, "Class distribution \\(classification\\):")
  expect_match(output, "Class distribution \\(validation\\):")
})

# Tests for print.train_models_rf_xgboost

test_that("print.train_models_rf_xgboost displays model info", {
  train_result <- list(
    cv_folds = 5,
    prop_train = 0.7,
    motifs_used = c("ATCG", "GCTA", "TACG"),
    actuals_test = factor(c("A", "B", "A", "B")),
    model_comparison = data.frame(
      model = c("rf", "xgboost"),
      test_accuracy = c(0.85, 0.88),
      validation_accuracy = c(0.82, 0.86)
    ),
    best_model_test = "xgboost",
    best_model_validation = "xgboost"
  )
  class(train_result) <- "train_models_rf_xgboost"

  output <- capture_print(train_result)

  expect_match(output, "Train Models RF/XGBoost Summary")
  expect_match(output, "CV folds: 5")
  expect_match(output, "Train proportion: 0.7")
  expect_match(output, "Motifs used: 3")
  expect_match(output, "Classes: 2")
  expect_match(output, "Performance comparison:")
  expect_match(output, "Best model \\(test\\): xgboost")
  expect_match(output, "Best model \\(validation\\): xgboost")
})

# Tests for print.kmer_analysis_result

test_that("print.kmer_analysis_result displays analysis info without validation", {
  analysis <- list(
    unique_cluster_motifs = data.frame(motif = c("ATCG", "GCTA")),
    unique_class_motifs = data.frame(motif = c("TACG")),
    cluster_frequency_ranking = data.frame(motif = rep("X", 10)),
    cluster_to_class = c(c1 = "A", c2 = "B"),
    class_frequency_matrix = data.frame(A = 1:5, B = 2:6, motif = letters[1:5]),
    has_external_validation = FALSE
  )
  class(analysis) <- "kmer_analysis_result"

  output <- capture_print(analysis)

  expect_match(output, "K-mer Analysis Results")
  expect_match(output, "TRAINING DATA:")
  expect_match(output, "Cluster-specific motifs:\\s+2")
  expect_match(output, "Class-specific motifs:\\s+1")
  expect_match(output, "Total motifs analyzed:\\s+10")
  expect_match(output, "Total clusters:\\s+2")
  expect_match(output, "Total classes:\\s+2")
  expect_false(grepl("VALIDATION DATA:", output))
})

test_that("print.kmer_analysis_result displays validation info when present", {
  analysis <- list(
    unique_cluster_motifs = data.frame(motif = c("ATCG")),
    unique_class_motifs = data.frame(motif = c("GCTA")),
    cluster_frequency_ranking = data.frame(motif = rep("X", 5)),
    cluster_to_class = c(c1 = "A"),
    class_frequency_matrix = data.frame(A = 1:5, motif = letters[1:5]),
    has_external_validation = TRUE,
    validation = list(
      unique_class_motifs = data.frame(motif = c("TACG", "CGAT")),
      n_sequences = 20,
      n_classes = 2
    )
  )
  class(analysis) <- "kmer_analysis_result"

  output <- capture_print(analysis)

  expect_match(output, "VALIDATION DATA:")
  expect_match(output, "Class-specific motifs:\\s+2")
  expect_match(output, "Total sequences:\\s+20")
  expect_match(output, "Total classes:\\s+2")
})

# Tests for print.kmers_in_seq_result

test_that("print.kmers_in_seq_result displays search results", {
  search_result <- data.frame(
    motif = rep(c("ATCG", "GCTA"), each = 5),
    sequence_name = paste0("seq", 1:10),
    class = rep(c("A", "B"), 5),
    position_start = seq(1, 100, length.out = 10),
    position_end = seq(4, 103, length.out = 10),
    sequence_length = rep(150, 10)
  )
  class(search_result) <- c("kmers_in_seq_result", "data.frame")
  attr(search_result, "n_train_sequences") <- 8
  attr(search_result, "has_validation") <- TRUE
  attr(search_result, "n_validation_sequences") <- 2
  attr(search_result, "n_motifs") <- 2
  attr(search_result, "n_occurrences") <- 10
  attr(search_result, "elapsed_time") <- 1.23

  output <- capture_print(search_result)

  expect_match(output, "K-mer Search Results")
  expect_match(output, "Training sequences: 8")
  expect_match(output, "Validation sequences: 2")
  expect_match(output, "Motifs searched: 2")
  expect_match(output, "Total occurrences: 10")
  expect_match(output, "Elapsed time: 1.23 s")
})

test_that("print.kmers_in_seq_result truncates long output", {
  search_result <- data.frame(
    motif = rep("ATCG", 50),
    sequence_name = paste0("seq", 1:50),
    class = rep("A", 50),
    position_start = 1:50,
    position_end = 5:54,
    sequence_length = rep(100, 50)
  )
  class(search_result) <- c("kmers_in_seq_result", "data.frame")
  attr(search_result, "n_train_sequences") <- 50
  attr(search_result, "has_validation") <- FALSE
  attr(search_result, "n_motifs") <- 1
  attr(search_result, "n_occurrences") <- 50
  attr(search_result, "elapsed_time") <- 0.5

  output <- capture_print(search_result)

  expect_match(output, "and 40 more rows")
})

# Tests for print.seq_classification

test_that("print.seq_classification displays complete pipeline results", {
  conf_mat <- list(
    table = matrix(c(8, 2, 1, 9), nrow = 2,
                   dimnames = list(Prediction = c("A", "B"),
                                   Reference = c("A", "B"))),
    overall = c(Accuracy = 0.85, Kappa = 0.70,
                AccuracyLower = 0.75, AccuracyUpper = 0.95)
  )

  classification <- list(
    parameters = list(
      k = 6,
      dist_method = "euclidean",
      hom_thresh = 0.8,
      seq_per_class = 10,
      min_size = 50,
      n_motifs = 20,
      prop_train = 0.7,
      cv_folds = 5,
      external_validation_fasta_dir = NULL
    ),
    best_model = "xgboost",
    validation_accuracy = 0.85,
    confusion_matrix = conf_mat,
    model_comparison = data.frame(
      model = c("rf", "xgboost"),
      test_accuracy = c(0.83, 0.85),
      validation_accuracy = c(0.80, 0.85)
    ),
    kmer_analysis = list(
      unique_cluster_motifs = data.frame(motif = letters[1:5]),
      unique_class_motifs = data.frame(motif = letters[1:3]),
      cluster_frequency_ranking = data.frame(motif = letters[1:15]),
      has_external_validation = FALSE
    ),
    processing_time = as.difftime(120, units = "secs"),
    timestamp = as.POSIXct("2024-01-01 12:00:00")
  )
  class(classification) <- "seq_classification"

  output <- capture_print(classification)

  expect_match(output, "SEQUENCE CLASSIFICATION RESULTS")
  expect_match(output, "PIPELINE PARAMETERS")
  expect_match(output, "K-mer size: 6")
  expect_match(output, "Distance method: euclidean")
  expect_match(output, "Validation type: Internal")
  expect_match(output, "BEST MODEL PERFORMANCE")
  expect_match(output, "Model: xgboost")
  expect_match(output, "Validation Accuracy: 0.85")
  expect_match(output, "CONFUSION MATRIX")
  expect_match(output, "MODEL COMPARISON")
  expect_match(output, "K-MER ANALYSIS SUMMARY")
  expect_match(output, "Cluster-specific motifs: 5")
  expect_match(output, "PROCESSING INFORMATION")
  expect_match(output, "Total time:")
})

test_that("print.seq_classification handles external validation", {
  conf_mat <- list(
    table = matrix(c(5, 0, 0, 5), nrow = 2),
    overall = c(Accuracy = 1.0, Kappa = 1.0,
                AccuracyLower = 0.95, AccuracyUpper = 1.0)
  )

  classification <- list(
    parameters = list(
      k = 4,
      dist_method = "manhattan",
      hom_thresh = 0.9,
      seq_per_class = 5,
      min_size = 100,
      n_motifs = 10,
      prop_train = 0.8,
      cv_folds = 10,
      external_validation_fasta_dir = "/path/to/validation"
    ),
    best_model = "rf",
    validation_accuracy = 0.95,
    confusion_matrix = conf_mat,
    model_comparison = data.frame(model = "rf", test_accuracy = 0.95,
                                  validation_accuracy = 0.95),
    kmer_analysis = list(
      unique_cluster_motifs = data.frame(motif = "A"),
      unique_class_motifs = data.frame(motif = "B"),
      cluster_frequency_ranking = data.frame(motif = "C"),
      has_external_validation = TRUE,
      validation = list(
        unique_class_motifs = data.frame(motif = c("D", "E")),
        class_frequency_ranking = data.frame(motif = letters[1:8])
      )
    ),
    processing_time = as.difftime(60, units = "secs"),
    timestamp = as.POSIXct("2024-01-01 12:00:00")
  )
  class(classification) <- "seq_classification"

  output <- capture_print(classification)

  expect_match(output, "Validation type: External")
  expect_match(output, "Validation class-specific motifs: 2")
  expect_match(output, "Validation motifs analyzed: 8")
})

# Test that all print methods return invisibly

test_that("all S3 print methods return invisibly", {
  # Test a sample of methods
  obj1 <- structure(list(kmers = data.frame(A = 1, CLASS = "X"),
                         metadata = data.frame(sequence_name = "s", length = 1, class = "X")),
                    class = c("kmer_data", "list"))
  expect_invisible(print(obj1))

  obj2 <- structure(list(n_elements = 1, cluster_summary = data.frame(),
                         n_unassigned = 0, min_size = 1, hom_thresh = 0.8),
                    class = c("cluster_dendrogram_result", "list"))
  expect_invisible(print(obj2))

  obj3 <- structure(list(A = "motif1"), class = "select_motifs",
                    n = 1, k = 1, m = 1, case = "test")
  expect_invisible(print(obj3))
})
