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

The `R` directory contains R code and supporting input files for analyses performed in the study.

The analyses are organized into separate subdirectories. Each subdirectory contains the corresponding R script, input files where applicable, and a README describing the analysis and how to run it.

Currently included analyses:

* `R/GSEA` – Gene Set Enrichment Analysis (GSEA) using DESeq2 results and MSigDB Human Hallmark gene sets.
* `R/Proteomics_Missing_Value_Imputation` – missing-value imputation of normalized proteomics data using the `missForest` R package.

Additional R analyses associated with the study may be provided in separate subdirectories.

## Repository structure

```text
melanopsin-paper/
├── python/
│   └── Python code for calcium signal analysis
│
├── R/
│   ├── GSEA/
│   └── Proteomics_Missing_Value_Imputation/
│
└── README.md
```

For analysis-specific details, required input files, and output descriptions, please refer to the README within each analysis directory.
