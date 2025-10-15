test_that("count_kmers returns correct counts for simple sequence", {
  # Create a simple test sequence
  seq <- Biostrings::DNAString("AAACCCGGGTTT")

  # Count 3-mers
  result <- count_kmers(seq, k = 3)

  # Check structure
  expect_type(result, "integer")
  expect_named(result)
  expect_equal(length(result), 4^3)  # All possible 3-mers with 4 letters

  # Check specific counts
  expect_equal(result["AAA"], c(AAA = 1))
  expect_equal(result["CCC"], c(CCC = 1))
  expect_equal(result["GGG"], c(GGG = 1))
  expect_equal(result["TTT"], c(TTT = 1))
})

test_that("count_kmers handles character input", {
  seq_char <- Biostrings::DNAString("ACGT")
  result <- count_kmers(seq_char, k = 2)

  expect_type(result, "integer")
  expect_equal(result["AC"], c(AC = 1))
  expect_equal(result["CG"], c(CG = 1))
  expect_equal(result["GT"], c(GT = 1))
})

test_that("count_kmers includes all possible kmers even with zero counts", {
  # Very short sequence won't have all kmers
  seq <- Biostrings::DNAString("AAA")
  result <- count_kmers(seq, k = 2)

  # Should still have all 16 possible 2-mers
  expect_equal(length(result), 16)

  # Some should be zero
  expect_equal(result["AA"], c(AA = 2))  # Present in sequence
  expect_equal(result["TT"], c(TT = 0))  # Not present in sequence
  expect_equal(result["GC"], c(GC = 0))  # Not present in sequence
})

test_that("count_kmers works with different k values", {
  seq <- Biostrings::DNAString("ACGTACGT")

  # Test k=1
  result_k1 <- count_kmers(seq, k = 1)
  expect_equal(length(result_k1), 4)
  expect_equal(result_k1["A"], c(A = 2))

  # Test k=2
  result_k2 <- count_kmers(seq, k = 2)
  expect_equal(length(result_k2), 16)

  # Test k=4
  result_k4 <- count_kmers(seq, k = 4)
  expect_equal(length(result_k4), 256)
})

test_that("count_kmers handles edge cases", {
  # Empty-ish sequence
  seq_short <- Biostrings::DNAString("A")
  result <- count_kmers(seq_short, k = 2)
  expect_equal(length(result), 16)
  expect_true(all(result == 0))  # No 2-mers in single letter

  # Sequence exactly k length
  seq_exact <- Biostrings::DNAString("ACG")
  result <- count_kmers(seq_exact, k = 3)
  expect_equal(result["ACG"], c(ACG = 1))
})
