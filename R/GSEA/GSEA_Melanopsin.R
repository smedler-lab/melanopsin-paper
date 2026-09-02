# Aim: Performing GSEA
# Written by: Parvaneh Nikpour (Smedler Lab)
# Initial steps before GSEA (preparation step)

# Checking the working directory (folder)
getwd()

# Changing the working directory. Set the working directory to the folder
# containing the input files, if needed.
# setwd("path/to/your/working/directory")

# Re-checking the working directory (folder)
getwd()

# Loading library for reading excel files
library(readxl)

# Load the DESeq2 results from the Excel file.
# DESeq2 results for both 1 h and 12 h are provided as input files.
# The code below shows the analysis for the 1 h dataset.
# To perform the same analysis for 12 h, replace the 1 h input file with the
# corresponding 12 h file and run the same workflow.

deseq_results_1 <- read_excel("DESeq2_dds_1h_diff.xlsx")

# For the 12 h analysis, use:
# deseq_results_1 <- read_excel("DESeq2_dds_12h_diff.xlsx")

# View here
View(deseq_results_1)

# Remove rows with missing gene symbols or log2 fold-change values
deseq_output_clean <- deseq_results_1[
  !is.na(deseq_results_1$gene_symbols) &
    !is.na(deseq_results_1$log2FoldChange),
]

# Check that gene symbols and log2 fold-change values contain no missing values
anyNA(deseq_output_clean$gene_symbols)
anyNA(deseq_output_clean$log2FoldChange)

# Load the necessary library
library(fgsea)

# Prepare data: Create a named vector of gene statistics using log2 fold change
# and rank genes from highest to lowest log2 fold change
gene_list <- with(deseq_output_clean, setNames(log2FoldChange, gene_symbols))
gene_list <- sort(gene_list, decreasing = TRUE)

# Read Hallmark gene sets from the provided GMT file (MSigDB, Human Hallmark
# gene sets, gene symbols; accessed 10 Jan 2025;
# https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H)

gene_sets <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt")

# As GSEA relies on a permutation-based approach to estimate the significance of
# enrichment scores, Each time you run the analysis, the random permutations may
# slightly vary unless a random seed is explicitly set.
# Solution: Set a random seed before running fgsea. This ensures reproducibility.
set.seed(42)  # Use any fixed number for the seed

# Run GSEA
gsea_results <- fgsea(pathways = gene_sets, stats = gene_list)

# Explore results
head(gsea_results)

View(gsea_results)

# This function displays the internal structure of the gsea_results object
str(gsea_results)

# Loads the data.table package into your R session, making its functions
# (like fwrite()) available.
library(data.table)

# Save the GSEA results as a CSV file
fwrite(gsea_results, "Melanopsin_Hela_1h_GSEA_results_Fixed_Seed.csv")

########################################################################
#################### Visualizing enrichment plots ######################
########################################################################

# Load ggplot2

library(ggplot2)

# Visualize the enrichment plot for a selected significant term (padj < 0.05).

# Replace "HALLMARK_PATHWAY_NAME" below with the name of the pathway of interest.

plotEnrichment(
  pathway = gene_sets[["HALLMARK_PATHWAY_NAME"]],
  stats = gene_list,
  gseaParam = 1,
  ticksSize = 0.5
) +
  labs(title = "HALLMARK_PATHWAY_NAME") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

########################################################################
##################### Saving plots in EPS format #######################
########################################################################

# Save the enrichment plot for a selected significant term in EPS format.

# Replace "HALLMARK_PATHWAY_NAME" with the name of the pathway of interest.

cairo_ps(
  file = "HALLMARK_PATHWAY_NAME.eps",
  width = 7,
  height = 5,
  family = "serif",
  onefile = FALSE
)

# Generate the enrichment plot
print(
  plotEnrichment(
    pathway = gene_sets[["HALLMARK_PATHWAY_NAME"]],
    stats = gene_list,
    gseaParam = 1,
    ticksSize = 0.5
  ) +
    labs(title = "HALLMARK_PATHWAY_NAME") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
)

# Close the EPS device
dev.off() 

########################################################################
######################## Saving R workspace ############################
########################################################################

# Save all objects in the current R session's global environment (workspace).

save.image(file = "Melanopsin_Hela_1h_GSEA_results.RData")

