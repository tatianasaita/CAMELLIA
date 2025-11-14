# tests/testthat/test-select_motifs.R
# ===== HELPER FUNCTIONS =====

.create_mock_motif_cluster <- function(n_kmers = 15, n_clusters = 5) {
  kmers <- replicate(n_kmers, {
    paste0(sample(c("A", "C", "G", "T"), 8, replace = TRUE), collapse = "")
  })

  kmers <- unique(kmers)
  while (length(kmers) < n_kmers) {
    new_kmer <- paste0(sample(c("A", "C", "G", "T"), 8, replace = TRUE),
                       collapse = "")
    if (!new_kmer %in% kmers) {
      kmers <- c(kmers, new_kmer)
    }
  }

  kmers <- unique(kmers)[1:n_kmers]

  cluster_data <- matrix(
    sample(0:10, n_kmers * n_clusters, replace = TRUE),
    nrow = n_kmers,
    ncol = n_clusters,
    dimnames = list(
      kmers,
      paste0("Cluster_", seq_len(n_clusters))
    )
  )

  return(as.data.frame(cluster_data))
}

.create_mock_cluster_result <- function(n_clusters = 5) {

  classes <- c("ClassA", "ClassB", "ClassC")

  # Distribuições pré-definidas (simples e previsível)
  class_distributions <- c(
    "ClassA(10), ClassB(5)",
    "ClassB(15), ClassC(3)",
    "ClassC(12), ClassA(4)",
    "ClassA(8), ClassB(9), ClassC(2)",
    "ClassB(14), ClassC(6), ClassA(1)",
    "ClassA(11), ClassB(7), ClassC(5)",
    "ClassC(13), ClassA(3), ClassB(8)",
    "ClassB(10), ClassC(4), ClassA(6)"
  )

  # Repete o padrão se necessário
  class_distributions <- rep(class_distributions,
                             length.out = n_clusters)

  cluster_summary <- data.frame(
    cluster_id = seq_len(n_clusters),
    dominant_class = rep(classes, length.out = n_clusters),
    class_distribution = class_distributions[1:n_clusters],
    stringsAsFactors = FALSE
  )

  list(cluster_summary = cluster_summary)
}
# ===== INPUT VALIDATION TESTS =====

test_that("Rejects non-numeric n", {
  motif_cluster <- .create_mock_motif_cluster()
  cluster_result <- .create_mock_cluster_result()

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = "5"),
    "must be a positive integer",
    ignore.case = TRUE
  )
})

test_that("Rejects non-positive n", {
  motif_cluster <- .create_mock_motif_cluster()
  cluster_result <- .create_mock_cluster_result()

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 0),
    "must be a positive integer",
    ignore.case = TRUE
  )
})

test_that("Rejects non-integer n", {
  motif_cluster <- .create_mock_motif_cluster()
  cluster_result <- .create_mock_cluster_result()

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 5.5),
    "must be a positive integer",
    ignore.case = TRUE
  )
})

test_that("Rejects non-data.frame motif_cluster", {
  cluster_result <- .create_mock_cluster_result()

  expect_error(
    select_motifs(list(), cluster_result, n = 5),
    "must be a data.frame",
    ignore.case = TRUE
  )
})

test_that("Rejects motif_cluster without rownames", {
  motif_cluster <- data.frame(Cluster_1 = c(1, 2, 3))
  cluster_result <- .create_mock_cluster_result()

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 5),
    "rownames",
    ignore.case = TRUE
  )
})

test_that("Rejects motif_cluster with only numeric rownames", {
  motif_cluster <- data.frame(Cluster_1 = c(1, 2, 3), row.names = c("1", "2", "3"))
  cluster_result <- .create_mock_cluster_result()

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 5),
    "meaningful rownames",
    ignore.case = TRUE
  )
})

test_that("Rejects motif_cluster with invalid column names", {
  motif_cluster <- data.frame(
    invalid_col = c(1, 2, 3),
    row.names = c("kmer1", "kmer2", "kmer3")
  )
  cluster_result <- .create_mock_cluster_result()

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 5),
    "Cluster_",
    ignore.case = TRUE
  )
})

test_that("Rejects non-list cluster_result", {
  motif_cluster <- .create_mock_motif_cluster()

  expect_error(
    select_motifs(motif_cluster, "not_a_list", n = 5),
    "must be a list",
    ignore.case = TRUE
  )
})


test_that("Rejects cluster_result without cluster_summary", {
  motif_cluster <- .create_mock_motif_cluster()
  cluster_result <- list(other_element = data.frame())

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 5),
    "cluster_summary",
    ignore.case = TRUE
  )
})

test_that("Rejects non-data.frame cluster_summary", {
  motif_cluster <- .create_mock_motif_cluster()
  cluster_result <- list(cluster_summary = list())

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 5),
    "must be a data.frame",
    ignore.case = TRUE
  )
})

test_that("Rejects cluster_summary missing required columns", {
  motif_cluster <- .create_mock_motif_cluster()
  cluster_result <- list(
    cluster_summary = data.frame(cluster_id = 1:3)
  )

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 5),
    "missing columns",
    ignore.case = TRUE
  )
})

test_that("Rejects non-logical verbose parameter", {
  motif_cluster <- .create_mock_motif_cluster()
  cluster_result <- .create_mock_cluster_result()

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 5, verbose = "yes"),
    "must be logical",
    ignore.case = TRUE
  )
})

# ===== CASE 1 TESTS (n < k) =====

test_that("CASE 1: Returns empty result when n < k", {
  n_clusters <- 5
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)
  )

  expect_s3_class(result, "empty_result")
  expect_length(result, 0)
})

test_that("CASE 1: Has correct attributes", {
  n_clusters <- 5
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)
  )

  expect_equal(attr(result, "case"), "CASE_1")
  expect_equal(attr(result, "n"), 2)
  expect_true(!is.null(attr(result, "reason")))
})

test_that("CASE 1: Returns invisible result", {
  n_clusters <- 5
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  expect_invisible(
    suppressMessages(
      select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)
    )
  )
})

# ===== CASE 2 TESTS (k <= n < m) =====

test_that("CASE 2: Triggered with k <= n < m", {
  n_clusters <- 8
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)
  )

  expect_equal(attr(result, "case"), "CASE_2")
})

test_that("CASE 2: Distributes motifs across classes", {
  n_clusters <- 8
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)
  )

  for (classe in names(result)) {
    expect_gt(length(result[[classe]]), 0)
  }
})

test_that("CASE 2: Total selected equals n", {
  n_clusters <- 8
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)
  )

  total_selected <- sum(sapply(result, length))
  expect_equal(total_selected, 5)
})

test_that("CASE 2: All motifs are valid k-mers", {
  n_clusters <- 8
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)
  )

  all_selected <- unlist(result)
  valid_kmers <- rownames(motif_cluster)

  expect_true(all(all_selected %in% valid_kmers))
})

test_that("CASE 2: No duplicate motifs", {
  n_clusters <- 8
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)
  )

  all_selected <- unlist(result)
  expect_equal(length(all_selected), length(unique(all_selected)))
})

test_that("CASE 2: Returns invisible result", {
  n_clusters <- 8
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  expect_invisible(
    suppressMessages(
      select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)
    )
  )
})

# ===== CASE 3 TESTS (n >= m) =====

test_that("CASE 3: Triggered with n >= m", {
  n_clusters <- 5
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 30, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 10, verbose = FALSE)
  )

  expect_equal(attr(result, "case"), "CASE_3")
})

test_that("CASE 3: Total selected equals n", {
  n_clusters <- 5
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 30, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 10, verbose = FALSE)
  )

  total_selected <- sum(sapply(result, length))
  expect_equal(total_selected, 10)
})

test_that("CASE 3: All motifs are valid k-mers", {
  n_clusters <- 5
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 30, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 10, verbose = FALSE)
  )

  all_selected <- unlist(result)
  valid_kmers <- rownames(motif_cluster)

  expect_true(all(all_selected %in% valid_kmers))
})

test_that("CASE 3: No duplicate motifs", {
  n_clusters <- 5
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 30, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 10, verbose = FALSE)
  )

  all_selected <- unlist(result)
  expect_equal(length(all_selected), length(unique(all_selected)))
})

test_that("CASE 3: Returns invisible result", {
  n_clusters <- 5
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 30, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  expect_invisible(
    suppressMessages(
      select_motifs(motif_cluster, cluster_result, n = 10, verbose = FALSE)
    )
  )
})

# ===== S3 METHOD TESTS =====

test_that("print.select_motifs works for CASE 2", {
  n_clusters <- 8
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)
  )

  expect_output(print(result), "Selected Motifs Summary")
})

test_that("print.select_motifs works for CASE 1", {
  n_clusters <- 5
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)
  )

  expect_output(print(result), "Empty")
})

test_that("summary.select_motifs works for CASE 2", {
  n_clusters <- 8
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)
  )

  expect_output(summary(result), "Detailed Summary")
})

test_that("print and summary return invisible", {
  n_clusters <- 8
  motif_cluster <- .create_mock_motif_cluster(n_kmers = 20, n_clusters = n_clusters)
  cluster_result <- .create_mock_cluster_result(n_clusters = n_clusters)

  suppressMessages(
    result <- select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)
  )

  expect_invisible(print(result))
  expect_invisible(summary(result))
})

# ===== HELPER FUNCTIONS =====

.create_mock_motif_cluster <- function(n_kmers = 15, n_clusters = 5) {
  kmers <- replicate(n_kmers, {
    paste0(sample(c("A", "C", "G", "T"), 8, replace = TRUE), collapse = "")
  })

  kmers <- unique(kmers)
  if (length(kmers) < n_kmers) {
    while (length(kmers) < n_kmers) {
      new_kmer <- paste0(sample(c("A", "C", "G", "T"), 8, replace = TRUE),
                         collapse = "")
      if (!new_kmer %in% kmers) {
        kmers <- c(kmers, new_kmer)
      }
    }
  }

  kmers <- unique(kmers)[1:n_kmers]

  cluster_data <- matrix(
    sample(0:10, n_kmers * n_clusters, replace = TRUE),
    nrow = n_kmers,
    ncol = n_clusters,
    dimnames = list(
      kmers,
      paste0("Cluster_", seq_len(n_clusters))
    )
  )

  return(as.data.frame(cluster_data))
}


