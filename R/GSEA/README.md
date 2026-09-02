# Gene Set Enrichment Analysis (GSEA)

**Analysis and R code:** Parvaneh Nikpour, Smedler Lab  
**Analysis date:** January 2025


This folder contains the R code and input files used to perform Gene Set Enrichment Analysis (GSEA) for the melanopsin HeLa cell analysis.

## Analysis

GSEA was performed using the `fgsea` R package. Differential expression results generated using DESeq2 were used as input for the analysis.

Genes were ranked from highest to lowest according to their log2 fold change (`log2FoldChange`). Gene Set Enrichment Analysis was then performed using the MSigDB Human Hallmark gene sets.

The analysis was performed separately for the 1 h and 12 h datasets using the same workflow.

## Files

* `GSEA_Melanopsin.R` – R script used to perform GSEA and generate enrichment plots.
* `5th_DESeq2_dds_1h_diff.xlsx` – DESeq2 differential expression results used as input for the 1 h analysis.
* `5th_DESeq2_dds_12h_diff.xlsx` – DESeq2 differential expression results used as input for the 12 h analysis.
* `h.all.v2024.1.Hs.symbols.gmt` – MSigDB Human Hallmark gene sets (gene symbols), version 2024.1.Hs, accessed 10 January 2025.

## Running the analysis

The R script demonstrates the analysis using the 1 h DESeq2 dataset:

```r
deseq_results_1 <- read_excel("5th_DESeq2_dds_1h_diff.xlsx")
```

To perform the corresponding 12 h analysis, replace the 1 h input file with:

```r
deseq_results_1 <- read_excel("5th_DESeq2_dds_12h_diff.xlsx")
```

The remaining analysis workflow is the same for both time points.

A fixed random seed (`set.seed(42)`) is used for reproducibility.

## Gene sets

Human Hallmark gene sets from the Molecular Signatures Database (MSigDB) were used:

`h.all.v2024.1.Hs.symbols.gmt`

MSigDB Human Hallmark gene sets:
https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H

## Output

The script generates GSEA results that can be exported as CSV files and includes code for visualizing selected significantly enriched pathways (`padj < 0.05`). Enrichment plots can also be exported in EPS format.
