# Aim: Differential abundance analysis of 1 h Proteomics data using limma
# Written by: Parvaneh Nikpour (Smedler Lab)
# Note: Proteomics data is available only for the 1h time point.

# Checking the working directory (folder)
getwd()

# Changing the working directory. Set the working directory to the folder
# containing the input files, if needed.
# setwd("path/to/your/working/directory")

# Re-checking the working directory (folder)
getwd()

# Load necessary libraries
library(readxl)
library(limma)
library(ggplot2)

########################################################################
################ Loading and formatting Proteomics data ################
########################################################################

# Load the normalized and imputed Proteomics data
Proteomics_imputed_data <- as.data.frame(
  read_excel(
    "Proteomics_imputed_data.xlsx",
    sheet = "Sheet1_transposed"
  )
)

# View the data
View(Proteomics_imputed_data)

# Check dimensions of the data
dim(Proteomics_imputed_data)

# Check whether there are any missing values
any(is.na(Proteomics_imputed_data))

# Use Protein_ID as row names
rownames(Proteomics_imputed_data) <-
  Proteomics_imputed_data$Protein_ID

# Remove the Protein_ID column after assigning it as row names
Proteomics_imputed_data <-
  Proteomics_imputed_data[, -1]

# Convert abundance values to numeric
Proteomics_imputed_data <- data.frame(
  lapply(
    Proteomics_imputed_data,
    function(x) as.numeric(as.character(x))
  ),
  check.names = FALSE,
  row.names = rownames(Proteomics_imputed_data)
)

# View the formatted Proteomics data
View(Proteomics_imputed_data)

# Check dimensions of the formatted data
dim(Proteomics_imputed_data)

# Apply log2 transformation
# Adding 1 prevents -Inf values if any abundance value is zero
exprs_log2 <- log2(Proteomics_imputed_data + 1)

# View the log2-transformed data
View(exprs_log2)

########################################################################
################ Loading sample metadata ###############################
########################################################################

# Load sample metadata
metadata <- as.data.frame(
  read_excel(
    "Experimental_design.xlsx",
    sheet = "Sample_metadata_renamed"
  )
)

# View sample metadata
View(metadata)

# Check dimensions of sample metadata
dim(metadata)

# View the first rows of sample metadata
head(metadata)

########################################################################
################ Preparing metadata for limma ##########################
########################################################################

# Reorder sample metadata to match the columns of the Proteomics matrix
metadata <- metadata[
  match(colnames(exprs_log2), metadata$Sample_IDs),
]

# Check that the sample order matches between the Proteomics data
# and sample metadata
all(colnames(exprs_log2) == metadata$Sample_IDs)

# Convert experimental variables to factors
metadata$Cell_type <- as.factor(metadata$Cell_type)
metadata$Stimulation <- as.factor(metadata$Stimulation)

# Set Melanopsin as the reference level for Cell_type
metadata$Cell_type <- relevel(
  metadata$Cell_type,
  ref = "Melanopsin"
)

# Set 1 h 5:120s as the reference level for Stimulation
metadata$Stimulation <- relevel(
  metadata$Stimulation,
  ref = "1 h 5:120s"
)

# Check factor levels
levels(metadata$Cell_type)
levels(metadata$Stimulation)

# View the prepared sample metadata
View(metadata)

########################################################################
################ Creating the limma design matrix ######################
########################################################################

# Create the design matrix including Cell_type and Stimulation
design <- model.matrix(
  ~ Cell_type + Stimulation,
  data = metadata
)

# View the design matrix
design

# Check coefficient names
colnames(design)

########################################################################
################ Fitting the limma model ###############################
########################################################################

# Fit a linear model to each protein
fit <- lmFit(
  exprs_log2,
  design
)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

########################################################################
################ Extracting differential abundance results #############
########################################################################

# Extract results for 1 h 5:60s versus 1 h 5:120s
# 1 h 5:120s is the reference level for Stimulation

results <- topTable(
  fit,
  coef = "Stimulation1 h 5:60s",
  number = Inf,
  adjust = "fdr"
)

# View the results
View(results)

# View the first rows of the results
head(results)

########################################################################
################ Saving differential abundance results #################
########################################################################

# Save the differential abundance results
write.csv(
  results,
  "Proteomics_1h_DE_Melanopsin_60s_vs_120s.csv"
)

########################################################################
################ Checking adjusted P values ############################
########################################################################

# Check whether there are any NA values in the adjusted P-value column
any_na <- any(is.na(results$adj.P.Val))

# Display the result
print(any_na)

########################################################################
################ Counting significant proteins ########################
########################################################################

# Count significantly differentially abundant proteins
# adjusted P value < 0.05
num_DEPs <- sum(
  results$adj.P.Val < 0.05,
  na.rm = TRUE
)

# Display the number of significant proteins
print(num_DEPs)

########################################################################
################ Preparing data for volcano plot #######################
########################################################################

# Create a copy of the limma results for plotting
dep_results <- results

# Rename columns for consistency with RNA-seq plotting code
dep_results$log2FoldChange <- dep_results$logFC
dep_results$padj <- dep_results$adj.P.Val

# View the prepared results
View(dep_results)

# Remove rows with NA adjusted P values
dep_results_noNA <- dep_results[
  !is.na(dep_results$padj),
]

# Create an empty label column
dep_results_noNA$label <- NA

# View the data used for the volcano plot
View(dep_results_noNA)

########################################################################
################ Counting up- and downregulated proteins ###############
########################################################################

# Identify significantly upregulated proteins
over_idx <- which(
  dep_results_noNA$padj < 0.05 &
    dep_results_noNA$log2FoldChange > 0
)

# Identify significantly downregulated proteins
under_idx <- which(
  dep_results_noNA$padj < 0.05 &
    dep_results_noNA$log2FoldChange < 0
)

# Count significantly upregulated proteins
over_count <- length(over_idx)

# Count significantly downregulated proteins
under_count <- length(under_idx)

# Display the counts
cat(
  "Number of overexpressed DEPs:",
  over_count,
  "\n"
)

cat(
  "Number of underexpressed DEPs:",
  under_count,
  "\n"
)

########################################################################
################ Volcano plot visualization ############################
########################################################################

# Define the adjusted P-value cutoff line
cutoff_lines <- list(
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "green"
  )
)

# Define colors for significant and non-significant proteins
overexpressed_color <- "red"
underexpressed_color <- "blue"
non_sig_color <- "grey"

# Create the volcano plot
volcano_plot <- ggplot(
  dep_results_noNA,
  aes(
    x = log2FoldChange,
    y = -log10(padj)
  )
) +
  geom_point(
    aes(
      color = factor(
        ifelse(
          padj < 0.05,
          ifelse(
            log2FoldChange > 0,
            "Upregulated",
            "Downregulated"
          ),
          "Non-Significant"
        ),
        levels = c(
          "Non-Significant",
          "Upregulated",
          "Downregulated"
        )
      )
    ),
    size = 2,
    alpha = 1
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "Upregulated" = overexpressed_color,
      "Downregulated" = underexpressed_color,
      "Non-Significant" = non_sig_color
    )
  ) +
  geom_text(
    aes(label = label),
    hjust = 1,
    vjust = 1,
    size = 3,
    check_overlap = TRUE,
    na.rm = TRUE,
    color = "black"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  ) +
  cutoff_lines +
  xlim(
    min(dep_results_noNA$log2FoldChange) - 0.5,
    max(dep_results_noNA$log2FoldChange) + 0.5
  ) +
  ylim(
    0,
    max(-log10(dep_results_noNA$padj)) + 1
  ) +
  labs(
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"Adjusted p-value")
  )

# Display the volcano plot
print(volcano_plot)

########################################################################
################ Saving the volcano plot in JPEG format ################
########################################################################

# Save the volcano plot as JPEG
ggsave(
  filename = "Proteomics_1h_volcano_plot.jpeg",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 600,
  units = "in"
)

########################################################################
################ Saving the volcano plot in EPS format #################
########################################################################

# Save the volcano plot as EPS
cairo_ps(
  file = "Proteomics_1h_volcano_plot.eps",
  width = 8,
  height = 6,
  onefile = FALSE
)

print(volcano_plot)

# Close the EPS device
dev.off()

########################################################################
######################## Saving R workspace ############################
########################################################################

# Save all objects in the current R session's global environment
save.image(file = "Proteomics_limma_DE.RData")
