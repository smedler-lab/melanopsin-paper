# Proteomics Missing-Value Imputation

**Analysis and R code:** [Parvaneh Nikpour](https://github.com/parvanehnikpour), Smedler Lab

This folder contains the R code and input file used for missing-value imputation of normalized proteomics data in the study:

**Nikpour P, Varas-Godoy M, Uhlén P, Smedler E.
Decoding calcium oscillation frequency in transcriptional regulation.
bioRxiv (2025).**

https://doi.org/10.1101/2025.10.10.676024

## Analysis

Missing values in the normalized proteomics data were imputed using the `missForest` R package.

The `Protein_Accession` column was used to identify proteins. The protein abundance data were converted to numeric values and transposed so that samples were represented as rows and proteins as columns before imputation.

Proteins with missing values across all samples were identified prior to imputation. Missing-value imputation was then performed using `missForest` with a maximum of 10 iterations (`maxiter = 10`) and 100 trees (`ntree = 100`).

No fixed random seed was used in the original imputation analysis. Because `missForest` uses a random forest-based procedure, the exact imputed values may therefore vary between independent runs.

## Files

* `Proteomics_Missing_Value_Imputation.R` – R script used for missing-value imputation of the normalized proteomics data.
* `Normalized_data.xlsx` – normalized proteomics data used as input for the imputation analysis.

## Running the analysis

The normalized proteomics data are read from the `Abundance_normalized` sheet of the input Excel file:

```r
xlsx_file <- "Normalized_data.xlsx"

Proteomics_normalized <- read_excel(
  xlsx_file,
  sheet = "Abundance_normalized"
)
```

The `Protein_Accession` column is used as the protein identifier. The abundance data are then transposed so that samples are represented as rows and proteins as columns.

Missing-value imputation is performed using `missForest`:

```r
Proteomics_missForest <- missForest(
  Proteomics_normalized_2_t,
  maxiter = 10,
  ntree = 100,
  variablewise = FALSE,
  decreasing = FALSE,
  verbose = FALSE,
  mtry = floor(sqrt(ncol(Proteomics_normalized_2_t))),
  replace = TRUE,
  parallelize = "no"
)
```

The number of missing values is checked before and after imputation to verify the imputation procedure.

## Output

The imputed proteomics data, including sample names, are exported as:

`Proteomics_imputed_data.xlsx`

The R workspace containing the objects generated during the analysis is saved as:

`Proteomics_Missing_Value_Imputation.RData`
