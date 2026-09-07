# Proteomics Principal Component Analysis (PCA)

**Analysis and R code:** [Parvaneh Nikpour](https://github.com/parvanehnikpour), Smedler Lab

This folder contains the R code and input files used for principal component analysis (PCA) of normalized and imputed proteomics data in the study:

**Nikpour P, Varas-Godoy M, Uhlén P, Smedler E.
Decoding calcium oscillation frequency in transcriptional regulation.
bioRxiv (2025).**

https://doi.org/10.1101/2025.10.10.676024

## Analysis

Principal component analysis (PCA) was performed on normalized and imputed proteomics data to examine variation among samples.

Proteomics data were available for the 1 h time point. Missing-value imputation was performed previously in a separate analysis using the `missForest` R package. The resulting imputed proteomics dataset was used directly as input for the PCA.

The input abundance matrix was formatted so that samples were represented as rows and proteins as columns. Protein abundance values were converted to numeric values before PCA.

PCA was performed in R using `prcomp()` without additional scaling:

```r
PCA_imputed <- prcomp(
  Proteomics_imputed_data,
  scale. = FALSE
)
```

The percentage of total variance explained by each principal component was calculated from the standard deviations of the principal components.

Both two-dimensional (2D) and three-dimensional (3D) PCA visualizations were generated. The 2D PCA plot displays PC1 and PC2, while the interactive 3D PCA plot displays PC1, PC2, and PC3.

Sample cell type was represented using different point shapes. Sample information, including cell type and stimulation condition, was included in the PCA data and is available interactively in the 3D plot.

## Files

* `Proteomics_PCA.R` – R script used for PCA and visualization of the normalized and imputed proteomics data.
* `Proteomics_imputed_data.xlsx` – normalized and imputed proteomics data used as input for PCA.
* `Experimental_design.xlsx` – sample metadata containing sample IDs, cell type, and stimulation information.

## Running the analysis

The normalized and imputed proteomics data are read directly from the input Excel file:

```r
Proteomics_imputed_data <- as.data.frame(
  read_excel(
    "Proteomics_imputed_data.xlsx"
  )
)
```

The `Sample` column is used to assign sample names as row names and is then removed from the abundance matrix. The remaining protein abundance values are converted to numeric values.

Sample metadata are read from the `Sample_metadata` sheet of `Experimental_design.xlsx`. The metadata include `Cell_type` and `Stimulation` information and are matched to the samples in the proteomics abundance matrix.

PCA is then performed on the imputed protein abundance matrix:

```r
PCA_imputed <- prcomp(
  Proteomics_imputed_data,
  scale. = FALSE
)
```

The percentage of variance explained by each principal component is calculated as:

```r
explained_variance <- PCA_imputed$sdev^2 /
  sum(PCA_imputed$sdev^2) * 100
```

A PCA data frame containing PC1, PC2, PC3, sample names, cell type, and stimulation information is generated for visualization.

The 2D PCA plot displays PC1 versus PC2 using `ggplot2`. Cell types are distinguished by point shape, and the percentages of variance explained by PC1 and PC2 are included in the axis labels.

The 3D PCA plot is generated using `plotly` and displays PC1, PC2, and PC3. The plot is interactive, allowing sample information to be viewed by hovering over individual data points.

## Output

The 2D PCA plot is exported in EPS format as:

`Proteomics_PCA_2D_plot.eps`

The interactive 3D PCA plot is exported as an HTML file:

`Proteomics_PCA_3D_plot.html`

The R workspace containing the objects generated during the analysis is saved as:

`Proteomics_PCA.RData`
