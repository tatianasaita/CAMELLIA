<!-- README.md is generated from README.Rmd. Please edit that file -->

# CAMELLIA

**C**ombining **A**lign**M**ent-free f**E**ature se**L**ection c**L**ass**I**fic**A**tion

CAMELLIA is an R package providing a complete pipeline for sequence classification using k-mer analysis, hierarchical clustering (or Affinity Propagation), motif selection, and machine learning models (Random Forest and XGBoost), including feature importance and SHAP-based model interpretation.

## Installation

CAMELLIA is not yet on CRAN. Install the development version directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("tatianasaita/CAMELLIA")
```

### Dependencies

Core functionality requires:

```r
install.packages(c("parallelDist", "randomForest"))
```

If you plan to use `ml_method = "xgb"` (including Step 9 — feature importance and SHAP analysis), the following packages are also required:

```r
install.packages(c("xgboost", "ggplot2", "shapviz"))
```

The pipeline checks for these dependencies automatically and will stop early with instructions if any are missing.

## Input Data

CAMELLIA expects FASTA files organized by class:

```
input_dir/
├── class1.fasta
├── class2.fasta
└── class3.fasta
```

Each FASTA file should contain sequences for one class:

```
>sequence1
ATCGATCGATCG
>sequence2
GCTAGCTAGCTA
```

## Complete Workflow

The full pipeline is run through a single function, `seq_classification()`, which executes an 8-step process (9 steps when `cluster_method = "dendrogram"` and `ml_method = "xgb"`, which adds a feature importance/SHAP step):

1. K-mer counting and data creation
2. Hierarchical clustering dendrogram (only if `cluster_method = "dendrogram"`)
3. Cluster validation and homogeneity assessment
4. Cluster motif calculation
5. Top discriminative motif selection
6. Train/validation/test dataset preparation
7. Model training (Random Forest and/or XGBoost)
8. K-mer analysis for cluster and class-specific motifs
9. Feature importance and SHAP analysis (only if `ml_method = "xgb"`)

### Basic usage — internal train/test split, dendrogram clustering

```r
result <- seq_classification(
  input_dir      = "path/to/fasta_files",
  k              = 6,
  ml_method      = "xgb",
  cluster_method = "dendrogram"
)
```

### Using Affinity Propagation clustering instead of a dendrogram

```r
result <- seq_classification(
  input_dir      = "path/to/fasta_files",
  k              = 6,
  seq_per_class  = 200,
  min_size       = 800,
  cluster_method = "apcluster",
  ap_r           = 2
)
```

### Using external FASTA files for training and/or testing

```r
result <- seq_classification(
  input_dir             = "CAMELLIA/extdata",
  k                     = 6,
  ml_method             = "xgb",
  cluster_method        = "dendrogram",
  external_train_fasta  = "CAMELLIA/datasets_train_test/train_01.fasta",
  external_test_fasta   = "CAMELLIA/datasets_train_test/test_01.fasta"
)
```

If only `external_test_fasta` is provided, the internal dataset is used for training and the external FASTA is used exclusively as the test set. If both `external_train_fasta` and `external_test_fasta` are provided, both splits come entirely from external files.

### Saving results to disk

By default, `seq_classification()` only prints progress and results to the console. To additionally save every plot generated during the pipeline (dendrogram, XGBoost feature importance, SHAP beeswarm, etc.) as individual PDF files, plus the full console log as a `.txt` file, provide an `output_dir`:

```r
result <- seq_classification(
  input_dir      = "path/to/fasta_files",
  k              = 6,
  ml_method      = "xgb",
  cluster_method = "dendrogram",
  output_dir     = "path/to/output"
)
```

This creates the following structure:

```
path/to/output/
└── Results/
    ├── plot_001.pdf
    ├── plot_002.pdf
    ├── plot_003.pdf
    ├── ...
    └── pipeline_log.txt
```

If `output_dir` is left as `NULL` (the default), no files are written and the behavior is unchanged.

## Key Parameters

| Parameter | Description | Default |
|----|----|----|
| `input_dir` | Path to directory containing FASTA files organized by class | — |
| `k` | K-mer size for sequence analysis | `6` |
| `dist_method` | Distance metric used for the dendrogram (see `parallelDist::parDist`) | `"euclidean"` |
| `hom_thresh` | Homogeneity threshold for cluster validation (0–1) | `0.8` |
| `seq_per_class` | Number of sequences to select per class for classification | `200` |
| `min_size` | Minimum sequence length | `800` |
| `n_motifs` | Number of top motifs to select per class | `256` |
| `prop_train` | Proportion of data used for training (0–1) | `0.7` |
| `cv_folds` | Number of cross-validation folds | `10` |
| `external_train_fasta` | Optional path to an external FASTA file for training | `NULL` |
| `external_test_fasta` | Optional path to an external FASTA file for testing | `NULL` |
| `kmer_analysis_threshold` | Threshold used in the k-mer analysis step | `0.8` |
| `cluster_method` | Clustering method: `"dendrogram"` or `"apcluster"` | `"dendrogram"` |
| `ap_r` | Distance metric parameter used by Affinity Propagation | `2` |
| `ml_method` | Classification model: `"xgb"` or `"rf"` | `"xgb"` |
| `verbose` | If `TRUE`, prints progress and intermediate results | `TRUE` |
| `output_dir` | Optional path where a `Results/` folder with PDF plots and a log file will be created | `NULL` |

## Output Object

`seq_classification()` returns a list of class `"seq_classification"` containing:

- `classification_results`: full model training results (Random Forest and/or XGBoost)
- `model_metrics`: performance metrics for the trained model(s)
- `predictions_train` / `predictions_test`: predictions on training and test sets
- `confusion_matrix_train` / `confusion_matrix_test`: confusion matrices
- `kmer_analysis`: k-mer analysis results (cluster and class-specific motifs)
- `feature_importance`: XGBoost feature importance table and SHAP results (only if `ml_method = "xgb"`)
- `processing_time`: total pipeline execution time
- `parameters`: list of all input parameters used
- `intermediate_results`: intermediate results from all pipeline steps (only if `verbose = TRUE`), including `data_result`, `dend_result`, `cluster_result`, `motif_result`, `select_motifs_result`, `datasets_traintest`, `classification_result`, `kmer_analysis`, and `feature_importance`
- `timestamp`: execution timestamp

## Performance Tips

1. **K-mer size**: Larger `k` values (5–7) capture more specific sequence patterns but increase computation time.
2. **Homogeneity threshold** (`hom_thresh`): Lower values (0.6–0.7) create more clusters; higher values (0.8–0.9) create fewer, purer clusters.
3. **Motif selection** (`n_motifs`): Start with a moderate number (e.g. 100–256) and increase if classification accuracy is low.
4. **Cross-validation** (`cv_folds`): Use 5–10 folds for more reliable performance estimates.
5. **ml_method**: Feature importance and SHAP analysis (Step 9) are only available for `ml_method = "xgb"`.

## Citation

[![DOI](https://zenodo.org/badge/1068178884.svg)](https://doi.org/10.5281/zenodo.21366714)



## Getting Help

For more details, see the package vignette:

```r
vignette("camellia-workflow", package = "CAMELLIA")
```

Or check the function documentation:

```r
?seq_classification
```

## License

GPL-3

## Issues

Report bugs and request features at:
<https://github.com/tatianasaita/CAMELLIA/issues>
