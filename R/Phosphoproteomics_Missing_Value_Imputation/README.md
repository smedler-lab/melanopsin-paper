# Phosphoproteomics Missing-Value Imputation

**Analysis and R code:** [Parvaneh Nikpour](https://github.com/parvanehnikpour), Smedler Lab

This folder contains the R code and input file used for missing-value imputation of normalized phosphoproteomics data in the study:

**Nikpour P, Varas-Godoy M, Uhlén P, Smedler E.
Decoding calcium oscillation frequency in transcriptional regulation.
bioRxiv (2025).**

https://doi.org/10.1101/2025.10.10.676024

## Analysis

Missing values in the normalized phosphoproteomics data were imputed using the `missForest` R package.

Phosphoproteomics data were available for the 1 h time point. Before this analysis, 49 phosphopeptide rows assigned to more than one master protein were removed from the dataset.

The `Modifications_in_Master_Proteins` column was used to identify phosphopeptides. Unique row names were generated using `make.unique()` to preserve phosphopeptides with duplicated modification identifiers while ensuring unique feature names.

The phosphopeptide abundance data were converted to numeric values and transposed so that samples were represented as rows and phosphopeptides as columns before imputation.

Samples and phosphopeptides with missing values across all entries were checked prior to imputation.

Missing-value imputation was performed using `missForest` with a maximum of 10 iterations (`maxiter = 10`) and 100 trees (`ntree = 100`). The number of variables randomly sampled at each split was defined as the square root of the number of phosphopeptide columns.

The mean abundance and standard deviation per sample were compared before and after imputation. The total number of missing values was also checked before and after imputation to verify the imputation procedure.

## Files

* `Phosphoproteomics_Missing_Value_Imputation.R` – R script used for missing-value imputation of the normalized phosphoproteomics data.
* `Phosphoproteomics_normalized_data.xlsx` – normalized and cleaned phosphoproteomics data used as input for the imputation analysis.

## Running the analysis

The normalized and cleaned phosphoproteomics data are read from the `Abundance_normalized` sheet of the input Excel file:

```r
Phosphoproteomics_normalized <- as.data.frame(
  read_excel(
    "Phosphoproteomics_normalized_data.xlsx",
    sheet = "Abundance_normalized"
  )
)
```

The `Modifications_in_Master_Proteins` column is used as the phosphopeptide identifier. Unique row names are generated using `make.unique()`:

```r
rownames(Phosphoproteomics_normalized) <- make.unique(
  as.character(
    Phosphoproteomics_normalized$Modifications_in_Master_Proteins
  )
)
```

The identifier column is then removed, the phosphopeptide abundance values are converted to numeric values, and the abundance matrix is transposed so that samples are represented as rows and phosphopeptides as columns.

Missing-value imputation is performed using `missForest`:

```r
Phosphoproteomics_missForest <- missForest(
  Phosphoproteomics_normalized_t,
  maxiter = 10,
  ntree = 100,
  variablewise = FALSE,
  decreasing = FALSE,
  verbose = FALSE,
  mtry = floor(
    sqrt(ncol(Phosphoproteomics_normalized_t))
  ),
  replace = TRUE,
  classwt = NULL,
  cutoff = NULL,
  strata = NULL,
  sampsize = NULL,
  nodesize = NULL,
  maxnodes = NULL,
  xtrue = NA,
  parallelize = c("no", "variables", "forests")
)
```

The imputed phosphoproteomics abundance matrix is extracted from the `missForest` output:

```r
Phosphoproteomics_imputed <-
  Phosphoproteomics_missForest$ximp
```

The numbers of missing values before and after imputation are compared, together with sample-level mean abundance and standard deviation, to assess the imputation results.

Before export, the original phosphopeptide identifiers are retained as column names:

```r
colnames(Phosphoproteomics_imputed) <-
  colnames(Phosphoproteomics_normalized_t)
```

Sample names are then added as the first column of the final imputed dataset.

## Output

The normalized and imputed phosphoproteomics data, including sample names and phosphopeptide identifiers, are exported as:

`Phosphoproteomics_imputed_data.xlsx`

The R workspace containing the objects generated during the analysis is saved as:

`Phosphoproteomics_Missing_Value_Imputation.RData`
