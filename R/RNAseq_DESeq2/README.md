# RNA-seq Differential Expression Analysis (DESeq2)

**Analysis and R code:** [Parvaneh Nikpour](https://github.com/parvanehnikpour), Smedler Lab

This folder contains the R code and input files used for differential expression analysis of RNA-seq data in the study:

**Nikpour P, Varas-Godoy M, Uhlén P, Smedler E.  
Decoding calcium oscillation frequency in transcriptional regulation.  
bioRxiv (2025).**

https://doi.org/10.1101/2025.10.10.676024

## Analysis

Differential expression analysis of the RNA-seq data was performed using the `DESeq2` R package.

The RNA-seq count matrix was prepared using gene symbols as row identifiers. Count values were converted to numeric values and rounded to integers before creating a `DESeqDataSet`.

Sample information was obtained from the accompanying metadata file. `Cell_type` and `Stimulation` were included in the DESeq2 design:

```r
design = ~ Cell_type + Stimulation
```

Differential expression analysis was performed separately for the 1 h and 12 h stimulation groups.

The R script provided in this folder demonstrates the complete workflow for the **1 h comparison**. The same workflow was applied to the 12 h samples by replacing the corresponding 1 h stimulation groups with the 12 h stimulation groups.

Genes with an adjusted p-value (`padj`) below 0.05 were considered statistically significant.

## Files

* `RNAseq_DESeq2.R` – R script demonstrating the DESeq2 differential expression analysis for the 1 h comparison.
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

### 1 h analysis

For the 1 h analysis, samples corresponding to the two 1 h stimulation conditions are selected:

```r
dds_1h <- dds[
  dds$Stimulation %in% c("1h5_60s", "1h5_120s"),
]
```

DESeq2 is then run on the 1 h dataset:

```r
dds_1h <- DESeq(dds_1h)
```

Differential expression results are obtained using the following contrast:

```r
dds_1h_diff <- data.frame(
  results(
    dds_1h,
    contrast = c(
      "Stimulation",
      "1h5_60s",
      "1h5_120s"
    )
  )
)
```

This contrast reports the differential expression for `1h5_60s` relative to `1h5_120s`.

### 12 h analysis

The same workflow was applied independently to the 12 h samples by replacing:

```text
1h5_60s   →   12h5_60s
1h5_120s  →   12h5_120s
```

The full analysis code is therefore shown only for the 1 h comparison to avoid duplication.

## Differentially expressed genes

The DESeq2 output includes the estimated log2 fold change (`log2FoldChange`), p-value (`pvalue`), and adjusted p-value (`padj`) for each gene.

Genes satisfying:

```r
padj < 0.05
```

were considered statistically significant.

Rows with missing adjusted p-values were excluded where required for downstream visualization and counting of significant genes.

## Visualization

The script includes code for visualizing the differential expression results using a volcano plot.

Genes are classified according to their adjusted p-value and direction of log2 fold change as:

* **Overexpressed** – `padj < 0.05` and `log2FoldChange > 0`
* **Underexpressed** – `padj < 0.05` and `log2FoldChange < 0`
* **Non-significant** – `padj >= 0.05`

The volcano plot can also be exported in EPS format.

## Output

The 1 h DESeq2 differential expression results are exported as:

`DESeq2_dds_1h_diff.xlsx`

The same workflow generates the corresponding DESeq2 results for the 12 h analysis.

The 1 h volcano plot can be exported as:

`volcano_plot_1h.eps`

The R workspace containing the objects generated during the analysis is saved as:

`RNAseq_DESeq2_1h.RData`
