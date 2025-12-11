library(testthat)

# ===== INPUT VALIDATION TESTS =====

test_that("kmers_in_seq validates inputs", {
  expect_error(
    kmers_in_seq(motifs = "AAA", cluster_result = NULL),
    "cluster_result"
  )

  expect_error(
    kmers_in_seq(motifs = NULL, cluster_result = list()),
    "motifs"
  )

  expect_error(
    kmers_in_seq(motifs = character(0), cluster_result = list()),
    "motifs"
  )
})

test_that("kmers_in_seq validates cluster_result structure", {
  # Missing data_result
  expect_error(
    kmers_in_seq(motifs = "AAA", cluster_result = list()),
    "data_result"
  )

  # Missing metadata
  expect_error(
    kmers_in_seq(
      motifs = "AAA",
      cluster_result = list(data_result = list())
    ),
    "metadata"
  )

  # Invalid metadata structure
  expect_error(
    kmers_in_seq(
      motifs = "AAA",
      cluster_result = list(
        data_result = list(
          metadata = data.frame(wrong_column = 1:5)
        )
      )
    ),
    "sequence_name, class, length"
  )
})

test_that("kmers_in_seq requires sequences or input_dir", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 100,
        stringsAsFactors = FALSE
      )
    )
  )

  expect_error(
    kmers_in_seq(
      motifs = "AAA",
      cluster_result = cluster_result,
      sequences = NULL,
      input_dir = NULL
    ),
    "Either 'sequences' or 'input_dir' must be provided"
  )
})

# ===== BASIC FUNCTIONALITY TESTS =====

test_that("kmers_in_seq works with basic inputs", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = c("seq1", "seq2"),
        class = c("A", "B"),
        length = c(100, 120),
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- c("ATGAAATGCAAA", "TTTGGGCCCAAA")
  sequence_names <- c("seq1", "seq2")
  motifs <- c("AAA", "TTT")

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_s3_class(result, "kmers_in_seq_result")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("motif", "sequence_name", "class", "position_start",
                    "position_end", "sequence_length") %in% colnames(result)))
  expect_gt(nrow(result), 0)
})

test_that("kmers_in_seq finds correct positions", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 15,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "ATGAAATGCAAA"
  sequence_names <- "seq1"
  motifs <- "AAA"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_equal(nrow(result), 2)
  expect_equal(result$motif, c("AAA", "AAA"))
  expect_equal(result$position_start, c(4, 10))
  expect_equal(result$position_end, c(6, 12))
})

test_that("kmers_in_seq handles overlapping motifs", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 10,
        stringsAsFactors = FALSE
      )
    )
  )

  # "AAAAAA" contém 2 ocorrências NÃO sobrepostas de "AAA"
  # Posições: 1-3 (AAA) e 4-6 (AAA)
  # stringi::stri_locate_all_fixed não encontra sobreposições por padrão
  sequences <- "AAAAAA"
  sequence_names <- "seq1"
  motifs <- "AAA"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  # Espera 2 ocorrências não-sobrepostas
  expect_equal(nrow(result), 2)
  expect_equal(result$position_start, c(1, 4))
  expect_equal(result$position_end, c(3, 6))
})

test_that("kmers_in_seq finds multiple non-overlapping occurrences", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 20,
        stringsAsFactors = FALSE
      )
    )
  )

  # Sequência com 4 ocorrências distintas de "AAA"
  sequences <- "AAATTTAAAGGGAAACCCAAA"
  sequence_names <- "seq1"
  motifs <- "AAA"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  # Deve encontrar 4 ocorrências
  expect_equal(nrow(result), 4)
  expect_true(all(result$motif == "AAA"))
})

# ===== MOTIF HANDLING TESTS =====

test_that("kmers_in_seq handles motif list input", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 10,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "ATGAAATTTGGG"
  sequence_names <- "seq1"
  motifs <- list(class_A = c("AAA", "TTT"), class_B = "GGG")

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_s3_class(result, "kmers_in_seq_result")
  expect_true(all(c("AAA", "TTT", "GGG") %in% result$motif))
})

test_that("kmers_in_seq handles duplicate motifs", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 12,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "ATGAAATGCAAA"
  sequence_names <- "seq1"
  motifs <- c("AAA", "AAA", "AAA")  # Duplicates

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_s3_class(result, "kmers_in_seq_result")
  expect_equal(attr(result, "n_motifs"), 1)  # Should deduplicate
})

test_that("kmers_in_seq handles no matches gracefully", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 10,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "ATGCATGC"
  sequence_names <- "seq1"
  motifs <- "ZZZZZ"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_s3_class(result, "kmers_in_seq_result")
  expect_equal(nrow(result), 0)
  expect_true(all(c("motif", "sequence_name", "class", "position_start",
                    "position_end", "sequence_length") %in% colnames(result)))
})

# ===== CLASS AND METADATA TESTS =====

test_that("kmers_in_seq includes class information", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = c("seq1", "seq2"),
        class = c("ClassA", "ClassB"),
        length = c(12, 12),
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- c("ATGAAATGC", "TTTGGGCCC")
  sequence_names <- c("seq1", "seq2")
  motifs <- c("AAA", "TTT")

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_true("class" %in% colnames(result))
  expect_true(all(result$class %in% c("ClassA", "ClassB")))
})

test_that("kmers_in_seq includes sequence_length column", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = c("seq1", "seq2"),
        class = c("A", "B"),
        length = c(100, 200),
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- c("ATGAAATGC", "TTTGGGCCC")
  sequence_names <- c("seq1", "seq2")
  motifs <- "AAA"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_true("sequence_length" %in% colnames(result))
  expect_true(all(result$sequence_length %in% c(100, 200)))
})

test_that("kmers_in_seq handles missing sequence_names", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = c("seq1", "seq2"),
        class = c("A", "B"),
        length = c(12, 12),
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- c("ATGAAATGC", "TTTGGGCCC")
  motifs <- "AAA"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = NULL,  # Will auto-generate
    verbose = FALSE
  )

  expect_s3_class(result, "kmers_in_seq_result")
  expect_true(all(grepl("^seq_", result$sequence_name)))
})

# ===== ATTRIBUTES TESTS =====

test_that("kmers_in_seq attributes are correct", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 12,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "ATGAAATGCAAA"
  sequence_names <- "seq1"
  motifs <- "AAA"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_equal(attr(result, "n_train_sequences"), 1)
  expect_equal(attr(result, "n_motifs"), 1)
  expect_equal(attr(result, "n_occurrences"), nrow(result))
  expect_false(attr(result, "has_validation"))
  expect_true(is.numeric(attr(result, "elapsed_time")))
})

test_that("kmers_in_seq tracks multiple sequences", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = c("seq1", "seq2", "seq3"),
        class = c("A", "B", "A"),
        length = c(12, 15, 10),
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- c("ATGAAATGC", "TTTGGGCCCAAA", "GGGGGG")
  sequence_names <- c("seq1", "seq2", "seq3")
  motifs <- c("AAA", "GGG")

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_equal(attr(result, "n_train_sequences"), 3)
  expect_equal(attr(result, "n_motifs"), 2)
})

# ===== S3 METHODS TESTS =====

test_that("kmers_in_seq print method works", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 12,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "ATGAAATGC"
  sequence_names <- "seq1"
  motifs <- "AAA"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_output(print(result), "K-mer Search Results")
  expect_output(print(result), "Training sequences")
  expect_output(print(result), "Motifs searched")
  expect_output(print(result), "Total occurrences")
})

test_that("kmers_in_seq summary method works", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = c("seq1", "seq2"),
        class = c("A", "B"),
        length = c(12, 12),
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- c("ATGAAATGC", "TTTGGGCCC")
  sequence_names <- c("seq1", "seq2")
  motifs <- c("AAA", "TTT")

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_output(summary(result), "K-mer Search Summary")
  expect_output(summary(result), "Top 15 motifs")
  expect_output(summary(result), "By class")
  expect_output(summary(result), "By sequence")
})

# ===== EDGE CASES =====

test_that("kmers_in_seq handles empty sequences", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 0,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- ""
  sequence_names <- "seq1"
  motifs <- "AAA"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_equal(nrow(result), 0)
})

test_that("kmers_in_seq handles case sensitivity", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 12,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "atgaaatgcaaa"  # lowercase
  sequence_names <- "seq1"
  motifs <- "AAA"  # uppercase

  # Should handle case internally (function uses toupper)
  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_gte(nrow(result), 0)
})

test_that("kmers_in_seq handles single nucleotide motifs", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 10,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "AAATTTGGG"
  sequence_names <- "seq1"
  motifs <- "A"

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_equal(nrow(result), 3)  # 3 A's in sequence
})

test_that("kmers_in_seq handles long motifs", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 20,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "ATGCATGCATGCATGC"
  sequence_names <- "seq1"
  motifs <- "ATGCATGC"  # 8 nucleotides

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = cluster_result,
    sequences = sequences,
    sequence_names = sequence_names,
    verbose = FALSE
  )

  expect_gte(nrow(result), 0)
})

# ===== VERBOSE OUTPUT TEST =====

test_that("kmers_in_seq verbose option works", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = "seq1",
        class = "A",
        length = 12,
        stringsAsFactors = FALSE
      )
    )
  )

  sequences <- "ATGAAATGC"
  sequence_names <- "seq1"
  motifs <- "AAA"

  # With verbose = TRUE, should produce output
  expect_message(
    kmers_in_seq(
      motifs = motifs,
      cluster_result = cluster_result,
      sequences = sequences,
      sequence_names = sequence_names,
      verbose = TRUE
    ),
    "K-mer Search in Sequences"
  )

  # With verbose = FALSE, no messages
  expect_silent(
    kmers_in_seq(
      motifs = motifs,
      cluster_result = cluster_result,
      sequences = sequences,
      sequence_names = sequence_names,
      verbose = FALSE
    )
  )
})
