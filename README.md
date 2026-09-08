# Melanopsin paper

This repository contains the analysis code and supporting files associated with the study:

**Nikpour P, Varas-Godoy M, Uhlén P, Smedler E.  
Decoding calcium oscillation frequency in transcriptional regulation.  
bioRxiv (2025).**

https://doi.org/10.1101/2025.10.10.676024

The repository includes both Python and R code used for data analysis and generation of figures presented in the study.

## Python code

The `python` directory contains the Python code used for calcium signal analysis.

To run the Python code, open the repository in a dev container in VS Code. More information about dev containers is available here:

https://code.visualstudio.com/docs/devcontainers/containers

Then install the required Python packages using:

```bash
pip install -r python/requirements.txt
```

The analysis can then be replicated by running:

`python/calcium_signal_analysis.ipynb`

## R code

The `R` directory contains R code and supporting input files for the RNA-seq, proteomics, and phosphoproteomics analyses performed in the study.

The analyses are organized into separate subdirectories. Each subdirectory contains the corresponding R script, input files where applicable, and a README describing the analysis and how to run it.

Currently included analyses:

### RNA-seq

* `R/RNAseq_PCA` – Principal Component Analysis (PCA) of normalized RNA-seq count data.
* `R/RNAseq_DESeq2` – differential expression analysis of RNA-seq data using DESeq2.
* `R/GSEA` – Gene Set Enrichment Analysis (GSEA) using DESeq2 results and MSigDB Human Hallmark gene sets.

### Proteomics

* `R/Proteomics_Missing_Value_Imputation` – missing-value imputation of normalized proteomics data using the `missForest` R package.
* `R/Proteomics_PCA` – Principal Component Analysis (PCA) of the imputed proteomics data.
* `R/Proteomics_limma_DE` – differential expression analysis of proteomics data using `limma`.

### Phosphoproteomics

* `R/Phosphoproteomics_Missing_Value_Imputation` – missing-value imputation of normalized phosphoproteomics data.
* `R/Phosphoproteomics_PCA` – Principal Component Analysis (PCA) of the imputed phosphoproteomics data.
* `R/Phosphoproteomics_limma_DE` – differential expression analysis of phosphoproteomics data using `limma`.

For analysis-specific information, including input files, analysis parameters, and outputs, please refer to the README within each analysis directory.

## Repository structure

```text
melanopsin-paper/
├── python/
│   └── Python code for calcium signal analysis
│
├── R/
│   ├── RNAseq_PCA/
│   ├── RNAseq_DESeq2/
│   ├── GSEA/
│   │
│   ├── Proteomics_Missing_Value_Imputation/
│   ├── Proteomics_PCA/
│   ├── Proteomics_limma_DE/
│   │
│   ├── Phosphoproteomics_Missing_Value_Imputation/
│   ├── Phosphoproteomics_PCA/
│   └── Phosphoproteomics_limma_DE/
│
└── README.md
```
**Analysis and R code:** [Parvaneh Nikpour](https://github.com/parvanehnikpour), Smedler Lab
