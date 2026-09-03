# Aim: Principal Component Analysis (PCA) of normalized RNA-seq counts
# Written by: Parvaneh Nikpour (Smedler Lab)
# Initial steps before PCA (preparation step)

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

# Round the count data to integer values
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

# View the prepared sample metadata
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
################ PCA using normalized counts ###########################
########################################################################

# Estimate size factors and extract normalized counts
dds_normalized <- estimateSizeFactors(dds)
normalized_counts <- counts(dds_normalized, normalized = TRUE)

# Transpose normalized counts so that samples are rows
normalized_counts_t <- t(normalized_counts)

# Perform PCA
PCA_normalized <- prcomp(normalized_counts_t)

# View PCA summary
summary(PCA_normalized)

# Create a data frame containing PCA coordinates and sample information
pca_df <- data.frame(
  PC1 = PCA_normalized$x[, 1],
  PC2 = PCA_normalized$x[, 2],
  PC3 = PCA_normalized$x[, 3],
  Sample = colData(dds)$Sample_IDs,
  Cell_type = colData(dds)$Cell_type,
  Stimulation = colData(dds)$Stimulation
)

# View PCA data
View(pca_df)

########################################################################
################ PCA visualization ####################################
########################################################################

# Calculate the percentage of variance explained by each principal component
explained_variance <- PCA_normalized$sdev^2 /
  sum(PCA_normalized$sdev^2) * 100

# Define axis labels showing variance explained
x_label <- paste0(
  "Principal Component 1 (",
  round(explained_variance[1], 1),
  "%)"
)

y_label <- paste0(
  "Principal Component 2 (",
  round(explained_variance[2], 1),
  "%)"
)

# Create the PCA plot
pca_plot <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, shape = Cell_type)
) +
  geom_point(size = 4, aes(fill = Cell_type), color = "black") +
  scale_shape_manual(values = c(21, 17)) +
  scale_fill_manual(values = c("white", NA)) +
  labs(
    x = x_label,
    y = y_label
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    plot.title = element_blank()
  )

# Display the PCA plot
print(pca_plot)

########################################################################
################ Saving the PCA plot in EPS format #####################
########################################################################

cairo_ps(
  file = "RNAseq_PCA_plot.eps",
  width = 6,
  height = 4,
  family = "serif",
  onefile = FALSE
)

print(pca_plot)

# Close the EPS device
dev.off()

########################################################################
######################## Saving R workspace ############################
########################################################################

# Save all objects in the current R session's global environment
save.image(file = "RNAseq_PCA.RData")