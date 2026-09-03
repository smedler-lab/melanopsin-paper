# Aim: Differential expression analysis of RNA-seq data using DESeq2
# Written by: Parvaneh Nikpour (Smedler Lab)
# Initial steps before DESeq2 analysis (preparation step)

# Checking the working directory (folder)
getwd()

# Changing the working directory. Set the working directory to the folder
# containing the input files, if needed.
# setwd("path/to/your/working/directory")

# Re-checking the working directory (folder)
getwd()

# Load necessary libraries
library(readxl)
library(DESeq2)
library(tibble)
library(writexl)
library(ggplot2)

########################################################################
################ Loading and formatting RNA-seq count data ############
########################################################################

# Load the RNA-seq count data
RNAseq1 <- as.data.frame(
  read_excel(
    "Quants_Bulk_Melanopsin.xlsx",
    sheet = "sheet1"
  )
)

# View the data
View(RNAseq1)

# Use gene symbols as row names
rownames(RNAseq1) <- RNAseq1$Gene_symbols

# Remove the Gene_symbols column after assigning it as row names
RNAseq1 <- RNAseq1[, -1]

# Convert count values to numeric
RNAseq1 <- data.frame(
  lapply(RNAseq1, function(x) as.numeric(as.character(x))),
  check.names = FALSE,
  row.names = rownames(RNAseq1)
)

# Round count data to integer values
RNAseq1 <- round(RNAseq1, digits = 0)

# View the formatted count data
View(RNAseq1)

########################################################################
################ Loading sample metadata ###############################
########################################################################

# Load sample metadata
Sample_info <- as.data.frame(
  read_excel(
    "Samples_metadata.xlsx",
    sheet = "All"
  )
)

# View sample metadata
View(Sample_info)

# Create a data frame containing sample names in the same order as the
# RNA-seq count matrix
group <- data.frame(colnames(RNAseq1))
colnames(group) <- "Sample_IDs"

# Merge sample names with sample metadata while preserving sample order
group2 <- merge(
  x = group,
  y = Sample_info,
  by = "Sample_IDs",
  sort = FALSE
)

# Convert experimental variables to factors
group2$Cell_type <- as.factor(group2$Cell_type)
group2$Stimulation <- as.factor(group2$Stimulation)

# View the prepared metadata
View(group2)

########################################################################
################ Creating the DESeq2 dataset ###########################
########################################################################

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = RNAseq1,
  colData = group2,
  design = ~ Cell_type + Stimulation
)

########################################################################
################ Differential expression analysis: 1 h ################
########################################################################

# The code below demonstrates the differential expression analysis for
# the 1 h comparison.
# The same workflow was applied to the corresponding 12 h samples by
# replacing the 1 h stimulation groups with "12h5_60s" and "12h5_120s".

# Subset the DESeq2 dataset for the 1 h stimulation samples
dds_1h <- dds[
  dds$Stimulation %in% c("1h5_60s", "1h5_120s"),
]

# Run DESeq2
dds_1h <- DESeq(dds_1h)

# Obtain differential expression results for 1h5_60s versus 1h5_120s
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

# Summarize the results
summary(dds_1h_diff)

# Add gene symbols as a column
dds_1h_diff <- rownames_to_column(
  dds_1h_diff,
  var = "Gene_symbols"
)

# View differential expression results
View(dds_1h_diff)

########################################################################
################ Checking DESeq2 results ###############################
########################################################################

# Check whether adjusted p-values contain missing values
any(is.na(dds_1h_diff$padj))

# Count genes with adjusted p-value < 0.05
num_DEGs_1h <- sum(
  dds_1h_diff$padj < 0.05,
  na.rm = TRUE
)

# Display the number of significant genes
cat(
  "Number of genes with padj < 0.05:",
  num_DEGs_1h,
  "\n"
)

# Remove rows with missing adjusted p-values
dds_1h_diff_noNA <- dds_1h_diff[
  !is.na(dds_1h_diff$padj),
]

# Check that adjusted p-values no longer contain missing values
any(is.na(dds_1h_diff_noNA$padj))

########################################################################
################ Saving DESeq2 results #################################
########################################################################

# Save the complete 1 h DESeq2 results as an Excel file
write_xlsx(
  dds_1h_diff,
  "DESeq2_dds_1h_diff.xlsx"
)

########################################################################
################ Volcano plot ##########################################
########################################################################

# Define gene categories according to adjusted p-value and log2 fold change
dds_1h_diff_noNA$category <- ifelse(
  dds_1h_diff_noNA$padj < 0.05 &
    dds_1h_diff_noNA$log2FoldChange > 0,
  "Overexpressed",
  ifelse(
    dds_1h_diff_noNA$padj < 0.05 &
      dds_1h_diff_noNA$log2FoldChange < 0,
    "Underexpressed",
    "Non-Significant"
  )
)

# Create the volcano plot
volcano_plot <- ggplot(
  dds_1h_diff_noNA,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = category
  )
) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(
    name = "Gene Status",
    values = c(
      "Overexpressed" = "red",
      "Underexpressed" = "blue",
      "Non-Significant" = "grey"
    )
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )

# Display the volcano plot
print(volcano_plot)

########################################################################
################ Saving volcano plot in EPS format #####################
########################################################################

cairo_ps(
  file = "volcano_plot_1h.eps",
  width = 7,
  height = 5,
  family = "serif",
  onefile = FALSE
)

print(volcano_plot)

# Close the EPS device
dev.off()

########################################################################
######################## Saving R workspace ############################
########################################################################

# Save all objects in the current R session's global environment
save.image(file = "RNAseq_DESeq2_1h.RData")