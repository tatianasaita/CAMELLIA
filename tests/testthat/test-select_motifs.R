# tests/testthat/test-select_motifs.R
#
# Essential tests for select_motifs function following CRAN standards
# These tests cover:
# - Input validation
# - Case 1: n < k (impossible)
# - Case 2: k <= n < m (distribute by class)
# - Case 3: n >= m (distribute by cluster)
# - S3 methods (print, summary)
# - Edge cases

test_that("select_motifs validates inputs correctly", {
  # Invalid n
  expect_error(
    select_motifs(data.frame(), list(), n = -1),
    "'n' must be a positive integer"
  )

  expect_error(
    select_motifs(data.frame(), list(), n = 2.5),
    "'n' must be a positive integer"
  )

  expect_error(
    select_motifs(data.frame(), list(), n = "5"),
    "'n' must be a positive integer"
  )

  # Invalid motif_cluster
  expect_error(
    select_motifs(list(), list(), n = 5),
    "'motif_cluster' must be a non-empty data.frame"
  )

  expect_error(
    select_motifs(data.frame(), list(), n = 5),
    "'motif_cluster' must be a non-empty data.frame"
  )

  # Invalid cluster_result
  expect_error(
    select_motifs(data.frame(x = 1), "not_a_list", n = 5),
    "'cluster_result' must contain 'cluster_summary' element"
  )

  expect_error(
    select_motifs(data.frame(x = 1), list(wrong = "element"), n = 5),
    "'cluster_result' must contain 'cluster_summary' element"
  )

  # Missing required columns
  motif_cluster <- data.frame(
    Cluster_1 = c(0.5, 0.8),
    row.names = c("AAA", "TTT")
  )

  cluster_result <- list(
    cluster_summary = data.frame(
      cluster_id = 1
      # missing dominant_class and class_distribution
    )
  )

  expect_error(
    select_motifs(motif_cluster, cluster_result, n = 1),
    "'cluster_summary' missing required columns"
  )
})

test_that("select_motifs handles CASE 1: n < k (impossible)", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.5, 0.8, 0.3),
    Cluster_2 = c(0.6, 0.7, 0.4),
    row.names = c("AAA", "TTT", "GGG")
  )

  cluster_summary <- data.frame(
    cluster_id = c(1, 2),
    dominant_class = c("ClassA", "ClassB"),
    class_distribution = c("ClassA(10)", "ClassB(10)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  # Try to select 1 motif when there are 2 classes
  result <- select_motifs(motif_cluster, cluster_result, n = 1, verbose = FALSE)

  # Check structure
  expect_s3_class(result, "select_motifs")
  expect_s3_class(result, "empty_result")
  expect_length(result, 0)

  # Check attributes
  expect_equal(attr(result, "n"), 1)
  expect_equal(attr(result, "k"), 2)
  expect_equal(attr(result, "case"), "CASE_1")
  expect_match(attr(result, "reason"), "n < k")
})

test_that("select_motifs handles CASE 2: k <= n < m", {
  # Create test data: 2 classes, 4 clusters, select 3 motifs
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.5, 0.3, 0.2, 0.1),
    Cluster_2 = c(0.8, 0.6, 0.4, 0.3, 0.2),
    Cluster_3 = c(0.7, 0.9, 0.5, 0.4, 0.3),
    Cluster_4 = c(0.6, 0.8, 0.7, 0.5, 0.4),
    row.names = c("AAA", "TTT", "GGG", "CCC", "ATC")
  )

  cluster_summary <- data.frame(
    cluster_id = c(1, 2, 3, 4),
    dominant_class = c("ClassA", "ClassA", "ClassB", "ClassB"),
    class_distribution = c("ClassA(10)", "ClassA(8)", "ClassB(12)", "ClassB(9)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 3, verbose = FALSE)

  # Check structure
  expect_s3_class(result, "select_motifs")
  expect_type(result, "list")

  # Check attributes
  expect_equal(attr(result, "n"), 3)
  expect_equal(attr(result, "k"), 2)
  expect_equal(attr(result, "m"), 4)
  expect_equal(attr(result, "case"), "CASE_2")

  # Check total motifs selected
  total_selected <- sum(sapply(result, length))
  expect_equal(total_selected, 3)

  # Check all classes are represented
  expect_true("ClassA" %in% names(result))
  expect_true("ClassB" %in% names(result))

  # Check motifs are unique
  all_motifs <- unlist(result, use.names = FALSE)
  expect_equal(length(all_motifs), length(unique(all_motifs)))
})

test_that("select_motifs handles CASE 2 with remainder", {
  # 2 classes, 4 clusters, select 5 motifs (base=2, remainder=1)
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05),
    Cluster_2 = c(0.8, 0.6, 0.4, 0.2, 0.15, 0.1),
    Cluster_3 = c(0.95, 0.8, 0.6, 0.4, 0.2, 0.1),
    Cluster_4 = c(0.85, 0.75, 0.55, 0.35, 0.25, 0.15),
    row.names = c("AAA", "TTT", "GGG", "CCC", "ATC", "GCT")
  )

  cluster_summary <- data.frame(
    cluster_id = c(1, 2, 3, 4),
    dominant_class = c("ClassA", "ClassA", "ClassB", "ClassB"),
    class_distribution = c("ClassA(15)", "ClassA(10)", "ClassB(8)", "ClassB(7)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 5, verbose = FALSE)

  # Check total
  total_selected <- sum(sapply(result, length))
  expect_equal(total_selected, 5)

  # ClassA (more present) should get the extra motif: 3 motifs
  # ClassB should get: 2 motifs
  expect_equal(length(result$ClassA), 3)
  expect_equal(length(result$ClassB), 2)
})

test_that("select_motifs handles CASE 3: n >= m", {
  # 2 classes, 3 clusters, select 6 motifs (2 per cluster)
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05, 0.02),
    Cluster_2 = c(0.8, 0.6, 0.4, 0.2, 0.15, 0.1, 0.05),
    Cluster_3 = c(0.95, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05),
    row.names = c("AAA", "TTT", "GGG", "CCC", "ATC", "GCT", "TAG")
  )

  cluster_summary <- data.frame(
    cluster_id = c(1, 2, 3),
    dominant_class = c("ClassA", "ClassA", "ClassB"),
    class_distribution = c("ClassA(10)", "ClassA(8)", "ClassB(12)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 6, verbose = FALSE)

  # Check structure
  expect_s3_class(result, "select_motifs")
  expect_equal(attr(result, "case"), "CASE_3")

  # Check total
  total_selected <- sum(sapply(result, length))
  expect_equal(total_selected, 6)

  # Check motifs are unique
  all_motifs <- unlist(result, use.names = FALSE)
  expect_equal(length(all_motifs), length(unique(all_motifs)))
})

test_that("select_motifs handles CASE 3 with remainder and clear class order", {
  # 2 classes, 3 clusters, select 7 motifs (base=2, remainder=1)
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2),
    Cluster_2 = c(0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15),
    Cluster_3 = c(0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25),
    row.names = c("AAA", "TTT", "GGG", "CCC", "ATC", "GCT", "TAG", "CAT")
  )

  # ClassX appears 5 times (more present)
  # ClassY appears 3 times (less present)
  cluster_summary <- data.frame(
    cluster_id = c(1, 2, 3),
    dominant_class = c("ClassY", "ClassX", "ClassX"),
    class_distribution = c("ClassY(10),ClassX(8)",
                           "ClassX(15),ClassY(5)",
                           "ClassX(18),ClassY(2)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 7, verbose = FALSE)

  # Check total
  total_selected <- sum(sapply(result, length))
  expect_equal(total_selected, 7)

  # Base: 2 per cluster (6 total)
  # Cluster 1 (ClassY): 2 → ClassY gets 2
  # Cluster 2 (ClassX): 2 → ClassX gets 2
  # Cluster 3 (ClassX): 2 → ClassX gets 4 total

  # Remainder: 1 goes to ClassX (more present: 5 vs 3)
  # Final: ClassY=2, ClassX=5
  expect_equal(length(result$ClassY), 2)
  expect_equal(length(result$ClassX), 5)

  # Verify most present class got more motifs
  expect_true(length(result$ClassX) > length(result$ClassY))
})

test_that("select_motifs prioritizes class with more cluster presence", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.8, 0.7, 0.6),
    Cluster_2 = c(0.85, 0.75, 0.65, 0.55),
    Cluster_3 = c(0.95, 0.85, 0.75, 0.65),
    Cluster_4 = c(0.92, 0.82, 0.72, 0.62),
    Cluster_5 = c(0.88, 0.78, 0.68, 0.58),
    row.names = c("AAA", "TTT", "GGG", "CCC")
  )

  # Beta: em 2 clusters (4, 5)
  # Alpha: em 3 clusters (1, 2, 3)
  cluster_summary <- data.frame(
    cluster_id = c(1, 2, 3, 4, 5),
    dominant_class = c("Alpha", "Alpha", "Alpha", "Beta", "Beta"),
    class_distribution = c("Alpha(10)",
                           "Alpha(15)",
                           "Alpha(18)",
                           "Beta(20)",
                           "Beta(15)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  # CASE 2: k=2, n=3, m=5 → 3 < 5
  result <- select_motifs(motif_cluster, cluster_result, n = 3, verbose = FALSE)

  # Cada classe aparece uma vez por cluster na distribution
  # unique() preserva: ["Alpha", "Beta"]
  # table() conta: Alpha=3 clusters, Beta=2 clusters
  # Ordem: Alpha (3), Beta (2)
  classes_order <- attr(result, "classes_order")
  expect_equal(classes_order[1], "Alpha")
  expect_equal(classes_order[2], "Beta")

  # Distribution: base=1, remainder=1
  # Alpha (mais presente) gets extra
  expect_equal(length(result$Alpha), 2)
  expect_equal(length(result$Beta), 1)
  expect_equal(sum(sapply(result, length)), 3)
})

test_that("select_motifs gives priority to class present in more clusters", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.8, 0.7, 0.6),
    Cluster_2 = c(0.85, 0.75, 0.65, 0.55),
    Cluster_3 = c(0.95, 0.85, 0.75, 0.65),
    Cluster_4 = c(0.92, 0.82, 0.72, 0.62),
    Cluster_5 = c(0.88, 0.78, 0.68, 0.58),
    row.names = c("AAA", "TTT", "GGG", "CCC")
  )

  # Major: 4 clusters, Minor: 1 cluster
  cluster_summary <- data.frame(
    cluster_id = c(1, 2, 3, 4, 5),
    dominant_class = c("Major", "Major", "Major", "Major", "Minor"),
    class_distribution = c("Major(10)", "Major(15)", "Major(12)", "Major(18)", "Minor(20)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 3, verbose = FALSE)

  expect_equal(attr(result, "classes_order")[1], "Major")
  expect_equal(length(result$Major), 2)
  expect_equal(length(result$Minor), 1)
})

test_that("select_motifs handles single cluster", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.7, 0.5),
    row.names = c("AAA", "TTT", "GGG")
  )

  cluster_summary <- data.frame(
    cluster_id = 1,
    dominant_class = "ClassA",
    class_distribution = "ClassA(20)"
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)

  expect_s3_class(result, "select_motifs")
  expect_equal(length(result$ClassA), 2)
  expect_equal(attr(result, "k"), 1)
  expect_equal(attr(result, "m"), 1)
})

test_that("select_motifs handles equal n and m", {
  # n = m = 3 (exactly one motif per cluster)
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.7, 0.5, 0.3),
    Cluster_2 = c(0.8, 0.6, 0.4, 0.2),
    Cluster_3 = c(0.95, 0.75, 0.55, 0.35),
    row.names = c("AAA", "TTT", "GGG", "CCC")
  )

  cluster_summary <- data.frame(
    cluster_id = c(1, 2, 3),
    dominant_class = c("ClassA", "ClassA", "ClassB"),
    class_distribution = c("ClassA(10)", "ClassA(8)", "ClassB(12)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 3, verbose = FALSE)

  expect_equal(attr(result, "case"), "CASE_3")
  expect_equal(sum(sapply(result, length)), 3)
})

test_that("print method works correctly", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.7),
    Cluster_2 = c(0.8, 0.6),
    row.names = c("AAA", "TTT")
  )

  cluster_summary <- data.frame(
    cluster_id = c(1, 2),
    dominant_class = c("ClassA", "ClassB"),
    class_distribution = c("ClassA(10)", "ClassB(8)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)

  expect_output(print(result), "Selected Motifs Summary")
  expect_output(print(result), "Total motifs")
  expect_output(print(result), "Classes")
  expect_output(print(result), "Clusters")
  expect_output(print(result), "Case")
  expect_output(print(result), "Motifs per class")
})

test_that("CASE 1: returns empty result when n < k", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.7),
    Cluster_2 = c(0.8, 0.6),
    row.names = c("AAA", "TTT")
  )

  # k = 2 dominant classes
  cluster_summary <- data.frame(
    cluster_id = c(1, 2),
    dominant_class = c("ClassA", "ClassB"),
    class_distribution = c("ClassA(10)", "ClassB(15)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  # n=1 < k=2
  result <- select_motifs(motif_cluster, cluster_result, n = 1, verbose = FALSE)

  expect_s3_class(result, "empty_result")
  expect_length(result, 0)
  expect_equal(attr(result, "case"), "CASE_1")
  expect_equal(attr(result, "k"), 2)
  expect_equal(attr(result, "n"), 1)
})

test_that("print method works for empty result (CASE 1)", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.8, 0.7),
    Cluster_2 = c(0.85, 0.75, 0.65),
    Cluster_3 = c(0.95, 0.85, 0.75),
    row.names = c("AAA", "TTT", "GGG")
  )

  # 3 dominant classes
  cluster_summary <- data.frame(
    cluster_id = c(1, 2, 3),
    dominant_class = c("Alpha", "Beta", "Gamma"),
    class_distribution = c("Alpha(10)", "Beta(15)", "Gamma(12)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  # n=2 < k=3 → CASE 1
  result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)

  expect_output(print(result), "Empty")
  expect_output(print(result), "n < k")
  expect_output(print(result), "insufficient motifs")
})

test_that("summary method works for empty result (CASE 1)", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.8),
    Cluster_2 = c(0.85, 0.75),
    row.names = c("AAA", "TTT")
  )

  cluster_summary <- data.frame(
    cluster_id = c(1, 2),
    dominant_class = c("X", "Y"),
    class_distribution = c("X(10)", "Y(15)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 1, verbose = FALSE)

  expect_output(summary(result), "Empty")
  expect_output(summary(result), "n < k")
})

test_that("summary method works correctly", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.7, 0.5),
    Cluster_2 = c(0.8, 0.6, 0.4),
    row.names = c("AAA", "TTT", "GGG")
  )

  cluster_summary <- data.frame(
    cluster_id = c(1, 2),
    dominant_class = c("ClassA", "ClassB"),
    class_distribution = c("ClassA(10)", "ClassB(8)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 2, verbose = FALSE)

  expect_output(summary(result), "Detailed Summary")
  expect_output(summary(result), "Selected motifs by class")
  expect_output(summary(result), "AAA|TTT|GGG")
})

test_that("verbose parameter works correctly", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.7),
    row.names = c("AAA", "TTT")
  )

  cluster_summary <- data.frame(
    cluster_id = 1,
    dominant_class = "ClassA",
    class_distribution = "ClassA(10)"
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  # verbose = TRUE should print messages
  expect_message(
    select_motifs(motif_cluster, cluster_result, n = 1, verbose = TRUE),
    "Classes"
  )

  # verbose = FALSE should not print messages
  expect_silent(
    select_motifs(motif_cluster, cluster_result, n = 1, verbose = FALSE)
  )
})

test_that("select_motifs works when only one class is dominant", {
  motif_cluster <- data.frame(
    Cluster_1 = c(0.9, 0.8, 0.7, 0.6),
    Cluster_2 = c(0.85, 0.75, 0.65, 0.55),
    Cluster_3 = c(0.95, 0.85, 0.75, 0.65),
    row.names = c("AAA", "TTT", "GGG", "CCC")
  )

  cluster_summary <- data.frame(
    cluster_id = c(1, 2, 3),
    dominant_class = c("Primary", "Primary", "Primary"),
    class_distribution = c("Primary(20),Secondary(5)",
                           "Primary(18),Secondary(7)",
                           "Primary(22),Secondary(3)")
  )

  cluster_result <- list(cluster_summary = cluster_summary)

  result <- select_motifs(motif_cluster, cluster_result, n = 3, verbose = FALSE)

  expect_equal(length(result$Primary), 3)
  expect_false("Secondary" %in% names(result))
  expect_equal(length(result), 1)
  expect_equal(names(result), "Primary")
})

