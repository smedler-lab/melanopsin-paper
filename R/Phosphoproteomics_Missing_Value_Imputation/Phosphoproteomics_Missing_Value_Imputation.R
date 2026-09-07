# Aim: Missing-value imputation of normalized Phosphoproteomics data
# Written by: Parvaneh Nikpour (Smedler Lab)
# Note: Phosphoproteomics data is available only for the 1 h time point.
# Before this analysis, 49 phosphopeptide rows assigned to more than one
# master protein were removed from the dataset.

# Checking the working directory (folder)
getwd()

# Changing the working directory. Set the working directory to the folder
# containing the input files, if needed.
# setwd("path/to/your/working/directory")

# Re-checking the working directory (folder)
getwd()

# Load necessary libraries
library(readxl)
library(missForest)
library(writexl)

########################################################################
######### Loading and formatting normalized Phosphoproteomics data #####
########################################################################

# Load the normalized and cleaned Phosphoproteomics data
Phosphoproteomics_normalized <- as.data.frame(
  read_excel(
    "Phosphoproteomics_normalized_data.xlsx",
    sheet = "Abundance_normalized"
  )
)

# View the data
View(Phosphoproteomics_normalized)

# View the first rows of the data
head(Phosphoproteomics_normalized)

# Check dimensions of the data
dim(Phosphoproteomics_normalized)

# Check whether there are any missing values
any(is.na(Phosphoproteomics_normalized))

# Count the total number of missing values
sum(is.na(Phosphoproteomics_normalized))

########################################################################
################ Formatting Phosphoproteomics data #####################
########################################################################

# Assign unique row names using Modifications_in_Master_Proteins
# make.unique() preserves duplicated phosphopeptide identifiers by adding
# unique suffixes where needed
rownames(Phosphoproteomics_normalized) <- make.unique(
  as.character(
    Phosphoproteomics_normalized$Modifications_in_Master_Proteins
  )
)

# View the data
View(Phosphoproteomics_normalized)

# Check dimensions of the data
dim(Phosphoproteomics_normalized)

# Remove the Modifications_in_Master_Proteins column after assigning
# it as row names
Phosphoproteomics_normalized_2 <-
  Phosphoproteomics_normalized[, -1]

# View the abundance data
View(Phosphoproteomics_normalized_2)

# Check dimensions of the abundance data
dim(Phosphoproteomics_normalized_2)

# Convert abundance values to numeric
Phosphoproteomics_normalized_2 <- data.frame(
  lapply(
    Phosphoproteomics_normalized_2,
    function(x) as.numeric(as.character(x))
  ),
  check.names = FALSE,
  row.names = rownames(Phosphoproteomics_normalized_2)
)

# View the formatted abundance data
View(Phosphoproteomics_normalized_2)

# Check dimensions of the formatted abundance data
dim(Phosphoproteomics_normalized_2)

########################################################################
############ Transposing Phosphoproteomics abundance data ##############
########################################################################

# Transpose the abundance matrix so that samples are rows and
# phosphopeptides are columns
Phosphoproteomics_normalized_t <- as.data.frame(
  t(Phosphoproteomics_normalized_2)
)

# View the transposed data
View(Phosphoproteomics_normalized_t)

# Check dimensions of the transposed data
dim(Phosphoproteomics_normalized_t)

# Check whether there are any missing values
any(is.na(Phosphoproteomics_normalized_t))

# Count the total number of missing values
sum(is.na(Phosphoproteomics_normalized_t))

########################################################################
################ Checking completely missing features ##################
########################################################################

# Identify samples with missing values across all phosphopeptides
which(
  rowSums(is.na(Phosphoproteomics_normalized_t)) ==
    ncol(Phosphoproteomics_normalized_t)
)

# Identify phosphopeptides with missing values across all samples
which(
  colSums(is.na(Phosphoproteomics_normalized_t)) ==
    nrow(Phosphoproteomics_normalized_t)
)

# Count samples with missing values across all phosphopeptides
sum(
  rowSums(is.na(Phosphoproteomics_normalized_t)) ==
    ncol(Phosphoproteomics_normalized_t)
)

# Count phosphopeptides with missing values across all samples
sum(
  colSums(is.na(Phosphoproteomics_normalized_t)) ==
    nrow(Phosphoproteomics_normalized_t)
)

########################################################################
################ Missing-value imputation ##############################
########################################################################

# Perform missing-value imputation using missForest
Phosphoproteomics_missForest <- missForest(
  Phosphoproteomics_normalized_t,
  maxiter = 10,
  ntree = 100,
  variablewise = FALSE,
  decreasing = FALSE,
  verbose = FALSE,
  mtry = floor(
    sqrt(ncol(Phosphoproteomics_normalized_t))
  ),
  replace = TRUE,
  classwt = NULL,
  cutoff = NULL,
  strata = NULL,
  sampsize = NULL,
  nodesize = NULL,
  maxnodes = NULL,
  xtrue = NA,
  parallelize = c("no", "variables", "forests")
)

# Display warnings generated during imputation, if any
warnings()

# View the missForest output
View(Phosphoproteomics_missForest)

# View the imputed abundance matrix
View(Phosphoproteomics_missForest$ximp)

# Check the out-of-bag imputation error
Phosphoproteomics_missForest$OOBerror

########################################################################
################ Checking the imputed data #############################
########################################################################

# Extract the imputed abundance matrix
Phosphoproteomics_imputed <-
  Phosphoproteomics_missForest$ximp

# Calculate mean abundance per sample before imputation
row_means_before <- rowMeans(
  Phosphoproteomics_normalized_t,
  na.rm = TRUE
)

# Calculate mean abundance per sample after imputation
row_means_after <- rowMeans(
  Phosphoproteomics_imputed
)

# View mean abundances before and after imputation
print(row_means_before)
print(row_means_after)

# Create a data frame comparing mean abundance before and after imputation
row_mean_comparison <- data.frame(
  Before_Imputation = row_means_before,
  After_Imputation = row_means_after
)

# View the comparison
print(row_mean_comparison)

# Compare mean abundance distributions before and after imputation
boxplot(
  row_means_before,
  row_means_after,
  names = c("Before", "After"),
  main = "Row-wise Mean (Per Sample)",
  ylab = "Mean Intensity"
)

# Calculate standard deviation per sample before imputation
row_sd_before <- apply(
  Phosphoproteomics_normalized_t,
  1,
  sd,
  na.rm = TRUE
)

# Calculate standard deviation per sample after imputation
row_sd_after <- apply(
  Phosphoproteomics_imputed,
  1,
  sd
)

# Compare standard deviation distributions before and after imputation
boxplot(
  row_sd_before,
  row_sd_after,
  names = c("Before", "After"),
  main = "Row-wise Standard Deviation",
  ylab = "Standard Deviation"
)

# Count missing values before imputation
total_na_before <- sum(
  is.na(Phosphoproteomics_normalized_t)
)

# Count missing values after imputation
total_na_after <- sum(
  is.na(Phosphoproteomics_imputed)
)

# Display the numbers of missing values before and after imputation
cat(
  "Missing values before:",
  total_na_before,
  "\n"
)

cat(
  "Missing values after:",
  total_na_after,
  "\n"
)

########################################################################
################ Saving the imputed data ###############################
########################################################################

# Ensure the original phosphopeptide identifiers are retained
# as column names after imputation
colnames(Phosphoproteomics_imputed) <-
  colnames(Phosphoproteomics_normalized_t)

# Add sample names as a new column
Phosphoproteomics_imputed_data <- cbind(
  Sample = rownames(Phosphoproteomics_imputed),
  Phosphoproteomics_imputed
)

# View the final imputed data
View(Phosphoproteomics_imputed_data)

# Check the first column names
head(colnames(Phosphoproteomics_imputed_data))

# Export the imputed Phosphoproteomics data
write_xlsx(
  Phosphoproteomics_imputed_data,
  "Phosphoproteomics_imputed_data.xlsx"
)

########################################################################
######################## Saving R workspace ############################
########################################################################

# Save all objects in the current R session's global environment
save.image(
  file = "Phosphoproteomics_Missing_Value_Imputation.RData"
)
