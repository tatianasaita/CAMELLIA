
#Funções auxiliares para executar seq_classification.R


setwd("E:/TATIANA/CAMELLIA-JUL26") #ajustar o caminho da pasta com as funções

source("create_data.R")
source("create_dendrogram.R")
source("cluster_dendrogram_ap1.R") 
source("calculate_cluster_motifs.R")
source("select_motifs.R")
source("select_train_test.R")
source("train_model_xgboost_rf.R")
source("kmer_analysis.R")
source("internal-functions.R")
source("methods.R")
source("seq_classification.R") 
source("feature_importance.R") 

################################# Ajustar o caminho das pastas ####################


result_01 <- seq_classification(
  input_dir            = "D:/datasets_train_semJ/train_01",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "D:/datasets_test_semJ/test_01.fasta",
  output_dir           = "E:/TATIANA/CAMELLIA-JUL26"
)

result_02 <- seq_classification(
  input_dir            = "datasets_train/train_02",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "datasets_test/test_02.fasta"
)

result_03 <- seq_classification(
  input_dir            = "datasets_train/train_03",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "datasets_test/test_03.fasta"
)

result_04 <- seq_classification(
  input_dir            = "datasets_train/train_04",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "datasets_test/test_04.fasta"
)

result_05 <- seq_classification(
  input_dir            = "datasets_train/train_05",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "datasets_test/test_05.fasta"
)

result_06 <- seq_classification(
  input_dir            = "datasets_train/train_06",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "datasets_test/test_06.fasta"
)

result_07 <- seq_classification(
  input_dir            = "datasets_train/train_07",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "datasets_test/test_07.fasta"
)

result_08 <- seq_classification(
  input_dir            = "datasets_train/train_08",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "datasets_test/test_08.fasta"
)

result_09 <- seq_classification(
  input_dir            = "datasets_train/train_09",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "datasets_test/test_09.fasta"
)

result_10 <- seq_classification(
  input_dir            = "datasets_train/train_10",
  k                    = 6,
  ml_method            = "xgb",
  cluster_method       = "dendrogram",
  external_test_fasta  = "datasets_test/test_10.fasta" 
)