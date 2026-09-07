# Phosphoproteomics Principal Component Analysis (PCA)

**Analysis and R code:** [Parvaneh Nikpour](https://github.com/parvanehnikpour), Smedler Lab

This folder contains the R code and input files used for principal component analysis (PCA) of normalized and imputed phosphoproteomics data in the study:

**Nikpour P, Varas-Godoy M, Uhlén P, Smedler E.
Decoding calcium oscillation frequency in transcriptional regulation.
bioRxiv (2025).**

https://doi.org/10.1101/2025.10.10.676024

## Analysis

Principal component analysis (PCA) was performed on normalized and imputed phosphoproteomics data to examine variation among samples.

Phosphoproteomics data were available for the 1 h time point. Missing-value imputation was performed previously in a separate analysis, and the resulting imputed phosphoproteomics dataset was used directly as input for PCA.

The input abundance matrix was formatted so that samples were represented as rows and phosphopeptides as columns. The `Sample` column was used to assign sample names as row names and was then removed from the abundance matrix.

Phosphopeptide abundance values were converted to numeric values before PCA.

Sample metadata included `Cell_type` and `Stimulation` information. The sample names from the phosphoproteomics abundance matrix were matched to the metadata, and the experimental variables were converted to factors.

PCA was performed in R using `prcomp()` without additional scaling:

```r
PCA_imputed <- prcomp(
  Phosphoproteomics_imputed_data,
  scale. = FALSE
)
```

The percentage of total variance explained by each principal component was calculated from the standard deviations of the principal components.

Both two-dimensional (2D) and three-dimensional (3D) PCA visualizations were generated. The 2D PCA plot displays PC1 and PC2, while the interactive 3D PCA plot displays PC1, PC2, and PC3.

Cell types were distinguished by point shape in the PCA plots. In the 2D plot, Empty cells are shown as open circles and Melanopsin cells as triangles. The 3D PCA plot uses the same black-and-white representation and includes sample, cell type, and stimulation information in the interactive hover text.

## Files

* `Phosphoproteomics_PCA.R` – R script used for PCA and visualization of the normalized and imputed phosphoproteomics data.
* `Phosphoproteomics_imputed_data.xlsx` – normalized and imputed phosphoproteomics data used as input for PCA.
* `Experimental_design.xlsx` – sample metadata containing sample IDs, cell type, and stimulation information.

## Running the analysis

The normalized and imputed phosphoproteomics data are read from the input Excel file:

```r
Phosphoproteomics_imputed_data <- as.data.frame(
  read_excel(
    "Phosphoproteomics_imputed_data.xlsx"
  )
)
```

The `Sample` column is used to assign sample names as row names and is then removed from the abundance matrix:

```r
rownames(Phosphoproteomics_imputed_data) <-
  Phosphoproteomics_imputed_data$Sample

Phosphoproteomics_imputed_data <-
  Phosphoproteomics_imputed_data[, -1]
```

The phosphopeptide abundance values are converted to numeric values before PCA.

Sample metadata are read from the `Sample_metadata` sheet of `Experimental_design.xlsx`. Sample names are matched to the phosphoproteomics abundance matrix, and `Cell_type` and `Stimulation` are converted to factors.

PCA is performed on the imputed phosphoproteomics abundance matrix:

```r
PCA_imputed <- prcomp(
  Phosphoproteomics_imputed_data,
  scale. = FALSE
)
```

The percentage of variance explained by each principal component is calculated as:

```r
explained_variance <- PCA_imputed$sdev^2 /
  sum(PCA_imputed$sdev^2) * 100
```

A PCA data frame containing PC1, PC2, PC3, sample names, cell type, and stimulation information is then generated for visualization.

The 2D PCA plot displays PC1 versus PC2 using `ggplot2`. The percentage of variance explained by each principal component is included in the axis labels.

The 3D PCA plot is generated using `plotly` and displays PC1, PC2, and PC3. The plot is interactive, allowing sample information, cell type, and stimulation condition to be viewed by hovering over individual points.

## Output

The 2D PCA plot is exported in EPS format as:

`Phosphoproteomics_PCA_2D_plot.eps`

The interactive 3D PCA plot is exported as an HTML file:

`Phosphoproteomics_PCA_3D_plot.html`

The R workspace containing the objects generated during the analysis is saved as:

`Phosphoproteomics_PCA.RData`
