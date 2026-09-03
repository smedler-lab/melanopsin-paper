# RNA-seq Principal Component Analysis (PCA)

**Analysis and R code:** [Parvaneh Nikpour](https://github.com/parvanehnikpour), Smedler Lab

This folder contains the R code and input files used for Principal Component Analysis (PCA) of RNA-seq data in the study:

**Nikpour P, Varas-Godoy M, Uhlén P, Smedler E.  
Decoding calcium oscillation frequency in transcriptional regulation.  
bioRxiv (2025).**

https://doi.org/10.1101/2025.10.10.676024

## Analysis

Principal Component Analysis (PCA) was performed on normalized RNA-seq count data to visualize the overall variation among samples.

The RNA-seq count matrix was prepared using gene symbols as row identifiers. Count values were converted to numeric values and rounded to integers before creating a `DESeqDataSet` using the `DESeq2` R package.

Sample information was obtained from the accompanying metadata file. `Cell_type` and `Stimulation` were included in the DESeq2 design.

Size-factor normalization was performed using `DESeq2`, and the normalized count matrix was transposed so that samples were represented as rows and genes as columns. PCA was then performed using the R function `prcomp()`.

## Files

* `RNAseq_PCA.R` – R script used for normalization, PCA, and visualization of the RNA-seq samples.
* `Quants_Bulk_Melanopsin.xlsx` – RNA-seq count data used as input for the analysis.
* `Samples_metadata.xlsx` – sample information used to define the experimental groups.

## Running the analysis

The RNA-seq count data are read from the `sheet1` sheet of the input Excel file:

```r
RNAseq1 <- as.data.frame(
  read_excel(
    "Quants_Bulk_Melanopsin.xlsx",
    sheet = "sheet1"
  )
)
```

Sample information is read from the `All` sheet of the metadata file:

```r
Sample_info <- as.data.frame(
  read_excel(
    "Samples_metadata.xlsx",
    sheet = "All"
  )
)
```

A DESeq2 dataset is then created using:

```r
dds <- DESeqDataSetFromMatrix(
  countData = RNAseq1,
  colData = group2,
  design = ~ Cell_type + Stimulation
)
```

Size factors are estimated and normalized counts are extracted:

```r
dds_normalized <- estimateSizeFactors(dds)

normalized_counts <- counts(
  dds_normalized,
  normalized = TRUE
)
```

The normalized count matrix is transposed, and PCA is performed using:

```r
normalized_counts_t <- t(normalized_counts)

PCA_normalized <- prcomp(
  normalized_counts_t,
  scale. = FALSE
)
```

The percentage of variance explained by each principal component is calculated from the PCA results and used in the axis labels of the PCA plot.

## Output

The script generates a PCA plot showing the distribution of RNA-seq samples according to the first two principal components.

The PCA plot can be exported in EPS format as:

`RNAseq_PCA_plot.eps`

The R workspace containing the objects generated during the analysis is saved as:

`RNAseq_PCA.RData`
