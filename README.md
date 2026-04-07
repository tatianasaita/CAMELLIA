<!-- README.md is generated from README.Rmd. Please edit that file -->

CAMELLIA expects FASTA files organized by class:

    input_dir/
    ├── class1.fasta
    ├── class2.fasta
    └── class3.fasta

Each FASTA file should contain sequences for one class:

    >sequence1
    ATCGATCGATCG
    >sequence2
    GCTAGCTAGCTA

## Complete Workflow

Basic usage with internal test split — apcluster method:

```r
result_seq_classification <- seq_classification_cent(
  input_dir     = "path/to/fasta_files",
  k             = 6,
  seq_per_class = 200,
  min_size      = 800,
  cluster_method = "apcluster",
  ap_r          = 2
)
```

## Key Parameters

| Parameter | Function | Description | Default |
|----|----|----|----|
| `k` | `create_data()` | K-mer size | \- |
| `dist_method` | `create_dendrogram()` | Distance metric | `"euclidean"` |
| `hom_thresh` | `cluster_dendrogram()` | Homogeneity threshold | `0.75` |
| `min_size` | `cluster_dendrogram()` | Minimum cluster size | `10` |
| `n` | `select_motifs()` | Number of motifs to select | `100` |
| `prop_train` | `train_models_rf_xgboost()` | Training proportion | `0.7` |
| `cv_folds` | `train_models_rf_xgboost()` | Cross-validation folds | `5` |

## Output Objects

All functions return S3 objects with specific classes:

- `camellia_data`: K-mer counts and metadata
- `camellia_dendrogram`: Hierarchical clustering results
- `camellia_clusters`: Cluster assignments and statistics
- `camellia_motifs`: Motif scores and rankings
- `camellia_models`: Trained classification models

Each object has `print()` and `summary()` methods for easy inspection.

## Performance Tips

1.  **K-mer size**: Larger k values (5-7) capture more specific patterns
    but increase computation time
2.  **Homogeneity threshold**: Lower values (0.6-0.7) create more
    clusters; higher values (0.8-0.9) create fewer, purer clusters
3.  **Motif selection**: Start with n=100 motifs; increase if
    classification accuracy is low
4.  **Cross-validation**: Use 5-10 folds for reliable performance
    estimates

## Citation

If you use CAMELLIA in your research, please cite:

    # Citation information here

## Getting Help

For more details, see the package vignette:

``` r
vignette("camellia-workflow", package = "CAMELLIA")
```

Or check the function documentation:

``` r
?create_data
?cluster_dendrogram
?train_models_rf_xgboost
```

## License

GPL-3

## Issues

Report bugs and request features at:
<https://github.com/tatianasaita/CAMELLIA/issues>
