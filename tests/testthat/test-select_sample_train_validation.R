
test_that("select_sample_train_validation rejects invalid cluster_result", {
  # Test 1: Missing 'class' column in metadata
  incomplete_metadata <- data.frame(
    sequence_name = c("seq1", "seq2"),
    length = c(400, 500)
    # Missing 'class' column
  )

  valid_kmers <- data.frame(
    kmer1 = c(1, 2),
    kmer2 = c(3, 4),
    row.names = c("seq1", "seq2")
  )

  cluster_result_invalid <- list(
    data_result = list(
      metadata = incomplete_metadata,
      kmers = valid_kmers
    )
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = cluster_result_invalid,
      k_per_class = 10,
      min_size = 300
    ),
    "Missing columns in metadata"
  )

  # Test 2: Missing 'data_result' key
  expect_error(
    select_sample_train_validation(
      cluster_result = list(wrong_key = "value"),
      k_per_class = 10,
      min_size = 300
    ),
    "must contain 'data_result'"
  )
})


test_that("select_sample_train_validation validates metadata and kmers", {
  invalid_result <- list(
    data_result = list(
      metadata = "not_a_dataframe",
      kmers = data.frame(a = 1)
    )
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = invalid_result,
      k_per_class = 10,
      min_size = 300
    ),
    "'metadata' must be a data.frame"
  )

  invalid_result2 <- list(
    data_result = list(
      metadata = data.frame(sequence_name = "seq1", length = 400, class = "A"),
      kmers = "not_a_dataframe"
    )
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = invalid_result2,
      k_per_class = 10,
      min_size = 300
    ),
    "'kmers' must be a data.frame"
  )
})


test_that("select_sample_train_validation validates k_per_class and min_size", {
  valid_metadata <- data.frame(
    sequence_name = c("seq1", "seq2"),
    length = c(400, 500),
    class = c("A", "B")
  )
  valid_kmers <- data.frame(a = 1:2, row.names = c("seq1", "seq2"))

  cluster_result <- list(
    data_result = list(
      metadata = valid_metadata,
      kmers = valid_kmers
    )
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = cluster_result,
      k_per_class = -5,
      min_size = 300
    ),
    "must be a positive integer"
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = cluster_result,
      k_per_class = 10.5,
      min_size = 300
    ),
    "must be a positive integer"
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = cluster_result,
      k_per_class = 10,
      min_size = -100
    ),
    "must be a positive integer"
  )
})


test_that("select_sample_train_validation requires necessary metadata columns", {
  incomplete_metadata <- data.frame(
    sequence_name = c("seq1", "seq2"),
    length = c(400, 500)
    # falta 'class'
  )
  valid_kmers <- data.frame(a = 1:2, row.names = c("seq1", "seq2"))

  cluster_result <- list(
    data_result = list(
      metadata = incomplete_metadata,
      kmers = valid_kmers
    )
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = cluster_result,
      k_per_class = 10,
      min_size = 300
    ),
    "Missing columns in metadata"
  )
})


test_that("select_sample_train_validation returns correct list structure", {
  # Criar dados válidos com múltiplas classes
  metadata <- data.frame(
    sequence_name = paste0("seq", 1:100),
    length = rep(c(400, 500, 600), length.out = 100),
    class = rep(c("A", "B", "C"), length.out = 100)
  )

  # Criar kmers com rownames correspondentes
  kmers <- data.frame(
    matrix(rnorm(100 * 50), nrow = 100),
    row.names = paste0("seq", 1:100)
  )

  cluster_result <- list(
    data_result = list(
      metadata = metadata,
      kmers = kmers
    )
  )

  result <- select_sample_train_validation(
    cluster_result = cluster_result,
    k_per_class = 10,
    min_size = 300
  )

  # Verificar estrutura do resultado
  expect_type(result, "list")
  expect_length(result, 4)
  expect_named(
    result,
    c(
      "true_labels_classification",
      "true_labels_validation",
      "classification_dataset",
      "validation_dataset"
    )
  )
})


test_that("select_sample_train_validation partitions data correctly", {
  metadata <- data.frame(
    sequence_name = paste0("seq", 1:50),
    length = rep(400, 50),
    class = rep(c("A", "B"), each = 25)
  )

  kmers <- data.frame(
    matrix(rnorm(50 * 30), nrow = 50),
    row.names = paste0("seq", 1:50)
  )

  cluster_result <- list(
    data_result = list(
      metadata = metadata,
      kmers = kmers
    )
  )

  result <- select_sample_train_validation(
    cluster_result = cluster_result,
    k_per_class = 5,
    min_size = 300
  )

  # Verificar contagem por classe no classification dataset
  classification_counts <- table(result$true_labels_classification$class)
  expect_equal(sum(classification_counts), 10) # 5 por classe * 2 classes

  # Verificar que não há sobreposição
  classification_names <- result$true_labels_classification$sequence_name
  validation_names <- result$true_labels_validation$sequence_name
  overlap <- intersect(classification_names, validation_names)
  expect_length(overlap, 0)

  # Verificar que as sequências correspondem aos dados de k-mers
  expect_true(
    all(classification_names %in% rownames(result$classification_dataset))
  )
})


test_that("select_sample_train_validation handles edge cases", {
  metadata <- data.frame(
    sequence_name = c("seq1", "seq2", "seq3"),
    length = c(400, 500, 600),
    class = c("A", "A", "B")
  )

  kmers <- data.frame(
    matrix(rnorm(3 * 10), nrow = 3),
    row.names = c("seq1", "seq2", "seq3")
  )

  cluster_result <- list(
    data_result = list(
      metadata = metadata,
      kmers = kmers
    )
  )

  # Redirecionar output temporariamente
  tmp <- tempfile()
  sink(tmp)

  result <- suppressWarnings(
    select_sample_train_validation(
      cluster_result = cluster_result,
      k_per_class = 100,
      min_size = 300
    )
  )

  sink() # Restaurar output
  file.remove(tmp) # Limpar arquivo temporário

  # Validações
  expect_type(result, "list")
  expect_true(nrow(result$true_labels_classification) > 0)
  expect_equal(nrow(result$true_labels_classification), 3)
})


test_that("select_sample_train_validation filters by minimum size", {
  metadata <- data.frame(
    sequence_name = paste0("seq", 1:10),
    length = c(100, 200, 300, 400, 500, 150, 350, 450, 550, 600),
    class = rep(c("A", "B"), each = 5)
  )

  kmers <- data.frame(
    matrix(rnorm(10 * 20), nrow = 10),
    row.names = paste0("seq", 1:10)
  )

  cluster_result <- list(
    data_result = list(
      metadata = metadata,
      kmers = kmers
    )
  )

  result <- select_sample_train_validation(
    cluster_result = cluster_result,
    k_per_class = 2,
    min_size = 300
  )

  # Verificar que todas as sequências selecionadas têm tamanho >= 300
  expect_true(
    all(result$true_labels_classification$length >= 300)
  )
})


test_that("select_sample_train_validation saves files", {
  metadata <- data.frame(
    sequence_name = paste0("seq", 1:20),
    length = rep(400, 20),
    class = rep(c("A", "B"), each = 10)
  )

  kmers <- data.frame(
    matrix(rnorm(20 * 15), nrow = 20),
    row.names = paste0("seq", 1:20)
  )

  cluster_result <- list(
    data_result = list(
      metadata = metadata,
      kmers = kmers
    )
  )

  # Usar temporary directory
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  temp_dir <- tempdir()
  setwd(temp_dir)

  result <- select_sample_train_validation(
    cluster_result = cluster_result,
    k_per_class = 3,
    min_size = 300
  )

  # Verificar se os arquivos foram criados
  expect_true(file.exists("true_labels_classification.RData"))
  expect_true(file.exists("true_labels_validation.RData"))
  expect_true(file.exists("classification_dataset.RData"))
  expect_true(file.exists("validation_dataset.RData"))
  expect_true(file.exists("backup_checkpoint1.RData"))
})

