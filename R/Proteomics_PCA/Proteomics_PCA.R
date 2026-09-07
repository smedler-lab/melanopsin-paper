# Aim: Principal Component Analysis (PCA) of normalized and imputed Proteomics data
# Written by: Parvaneh Nikpour (Smedler Lab)
# Note: Proteomics data is available only for the 1 h time point.
# Missing-value imputation was performed previously in a separate analysis.

# Checking the working directory (folder)
getwd()

# Changing the working directory. Set the working directory to the folder
# containing the input files, if needed.
# setwd("path/to/your/working/directory")

# Re-checking the working directory (folder)
getwd()

# Load necessary libraries
library(readxl)
library(ggplot2)
library(plotly)
library(htmlwidgets)

########################################################################
############ Loading and formatting imputed Proteomics data ############
########################################################################

# Load the normalized and imputed Proteomics data
Proteomics_imputed_data <- as.data.frame(
  read_excel(
    "Proteomics_imputed_data.xlsx"
  )
)

# View the data
View(Proteomics_imputed_data)

# Check dimensions of the data
dim(Proteomics_imputed_data)

# Check whether there are any missing values
any(is.na(Proteomics_imputed_data))

# Use sample names as row names
rownames(Proteomics_imputed_data) <-
  Proteomics_imputed_data$Sample

# Remove the Sample column after assigning it as row names
Proteomics_imputed_data <-
  Proteomics_imputed_data[, -1]

# Convert protein abundance values to numeric
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

########################################################################
################ Loading sample metadata ###############################
########################################################################

# Load sample metadata
Sample_info <- as.data.frame(
  read_excel(
    "Experimental_design.xlsx",
    sheet = "Sample_metadata"
  )
)

# View sample metadata
View(Sample_info)

# Create a data frame containing sample names in the same order as the
# Proteomics abundance matrix
group <- data.frame(rownames(Proteomics_imputed_data))
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
################ PCA using imputed Proteomics data #####################
########################################################################

# Perform PCA
# Samples are rows and proteins are columns
PCA_imputed <- prcomp(
  Proteomics_imputed_data,
  scale. = FALSE
)

# View PCA summary
summary(PCA_imputed)

# Calculate the percentage of variance explained by each principal component
explained_variance <- PCA_imputed$sdev^2 /
  sum(PCA_imputed$sdev^2) * 100

# Create a data frame containing PCA coordinates and sample information
pca_df <- data.frame(
  PC1 = PCA_imputed$x[, 1],
  PC2 = PCA_imputed$x[, 2],
  PC3 = PCA_imputed$x[, 3],
  Sample = rownames(PCA_imputed$x),
  Cell_type = group2$Cell_type,
  Stimulation = group2$Stimulation
)

# View PCA data
View(pca_df)

########################################################################
################ 2D PCA visualization #################################
########################################################################

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

# Create the 2D PCA plot
pca_plot_2d <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, shape = Cell_type)
) +
  geom_point(
    size = 4,
    aes(fill = Cell_type),
    color = "black"
  ) +
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
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    plot.title = element_blank()
  )

# Display the 2D PCA plot
print(pca_plot_2d)

########################################################################
################ Saving the 2D PCA plot in EPS format ##################
########################################################################

cairo_ps(
  file = "Proteomics_PCA_2D_plot.eps",
  width = 6,
  height = 4,
  family = "serif",
  onefile = FALSE
)

print(pca_plot_2d)

# Close the EPS device
dev.off()

########################################################################
################ 3D PCA visualization #################################
########################################################################

# Create an interactive 3D PCA plot using PC1, PC2, and PC3
pca_plot_3d <- plot_ly(
  pca_df,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  type = "scatter3d",
  mode = "markers",
  symbol = ~Cell_type,
  symbols = c(
    "Empty" = "circle-open",
    "Melanopsin" = "triangle-up"
  ),
  marker = list(
    size = 5,
    color = "black"
  ),
  text = ~paste0(
    "Sample: ", Sample,
    "<br>Cell type: ", Cell_type,
    "<br>Stimulation: ", Stimulation
  ),
  hoverinfo = "text"
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = paste0(
          "PC1 (",
          round(explained_variance[1], 1),
          "%)"
        ),
        showticklabels = FALSE
      ),
      yaxis = list(
        title = paste0(
          "PC2 (",
          round(explained_variance[2], 1),
          "%)"
        ),
        showticklabels = FALSE
      ),
      zaxis = list(
        title = paste0(
          "PC3 (",
          round(explained_variance[3], 1),
          "%)"
        ),
        showticklabels = FALSE
      )
    ),
    showlegend = TRUE
  )

# Display the 3D PCA plot
pca_plot_3d

########################################################################
################ Saving the 3D PCA plot as HTML ########################
########################################################################

# Save the interactive 3D PCA plot
saveWidget(
  pca_plot_3d,
  "Proteomics_PCA_3D_plot.html",
  selfcontained = TRUE
)

########################################################################
######################## Saving R workspace ############################
########################################################################

# Save all objects in the current R session's global environment
save.image(file = "Proteomics_PCA.RData")
