# Proteomics Differential Abundance Analysis

**Analysis and R code:** [Parvaneh Nikpour](https://github.com/parvanehnikpour), Smedler Lab

This folder contains the R code and input files used for differential abundance analysis of normalized and imputed proteomics data in the study:

**Nikpour P, Varas-Godoy M, Uhlén P, Smedler E.
Decoding calcium oscillation frequency in transcriptional regulation.
bioRxiv (2025).**

https://doi.org/10.1101/2025.10.10.676024

## Analysis

Differential protein abundance was analyzed using the `limma` R package.

Proteomics data were available for the 1 h time point. The normalized and imputed protein abundance data were log2-transformed after adding 1 to avoid infinite values for zero abundances.

Sample metadata included `Cell_type` and `Stimulation` as experimental variables. The metadata were reordered to match the sample order in the proteomics abundance matrix before differential abundance analysis.

`Cell_type` and `Stimulation` were treated as factors. `Melanopsin` was set as the reference level for `Cell_type`, and `1 h 5:120s` was set as the reference level for `Stimulation`.

The limma design matrix was defined as:

```r
design <- model.matrix(
  ~ Cell_type + Stimulation,
  data = metadata
)
```

A linear model was fitted to each protein using `lmFit()`, followed by empirical Bayes moderation using `eBayes()`.

Differential abundance was evaluated for `1 h 5:60s` relative to the reference stimulation condition `1 h 5:120s` using:

```r
results <- topTable(
  fit,
  coef = "Stimulation1 h 5:60s",
  number = Inf,
  adjust = "fdr"
)
```

P values were corrected for multiple testing using false discovery rate (FDR) adjustment. Proteins with an adjusted P value below 0.05 were considered significantly differentially abundant.

## Files

* `Proteomics_limma_DE.R` – R script used for differential abundance analysis of the proteomics data.
* `Proteomics_imputed_data.xlsx` – normalized and imputed proteomics data used as input for the limma analysis.
* `Experimental_design.xlsx` – sample metadata containing sample IDs, cell type, and stimulation information.

## Running the analysis

The normalized and imputed proteomics data are read from the `Sheet1_transposed` sheet of the input Excel file:

```r
Proteomics_imputed_data <- as.data.frame(
  read_excel(
    "Proteomics_imputed_data.xlsx",
    sheet = "Sheet1_transposed"
  )
)
```

The `Protein_ID` column is used as the protein identifier and assigned as row names. The remaining protein abundance values are converted to numeric values and log2-transformed:

```r
rownames(Proteomics_imputed_data) <-
  Proteomics_imputed_data$Protein_ID

Proteomics_imputed_data <-
  Proteomics_imputed_data[, -1]

exprs_log2 <- log2(Proteomics_imputed_data + 1)
```

Sample metadata are read from the `Sample_metadata_renamed` sheet of `Experimental_design.xlsx` and reordered to match the samples in the proteomics abundance matrix.

The experimental variables are converted to factors, and the reference levels are defined as:

```r
metadata$Cell_type <- relevel(
  metadata$Cell_type,
  ref = "Melanopsin"
)

metadata$Stimulation <- relevel(
  metadata$Stimulation,
  ref = "1 h 5:120s"
)
```

The limma model is then fitted using:

```r
design <- model.matrix(
  ~ Cell_type + Stimulation,
  data = metadata
)

fit <- lmFit(
  exprs_log2,
  design
)

fit <- eBayes(fit)
```

Differential abundance results are extracted for `1 h 5:60s` relative to `1 h 5:120s`. The number of significantly differentially abundant proteins is determined using an adjusted P-value threshold of 0.05.

A volcano plot is generated using the log2 fold change and adjusted P value. Significantly differentially abundant proteins with positive log2 fold changes are shown as upregulated, while those with negative log2 fold changes are shown as downregulated. No additional fold-change cutoff is applied.

## Output

The differential abundance results are exported as:

`Proteomics_1h_DE_Melanopsin_60s_vs_120s.csv`

The volcano plot is exported in JPEG and EPS formats as:

`Proteomics_1h_volcano_plot.jpeg`

`Proteomics_1h_volcano_plot.eps`

The R workspace containing the objects generated during the analysis is saved as:

`Proteomics_limma_DE.RData`
