# Aim: Imputation of missing values in proteomics data using missForest
# Written by: Parvaneh Nikpour (Smedler Lab)
# Initial steps before imputation (preparation step)

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

# Load the normalized proteomics data
xlsx_file <- "Normalized_data.xlsx"

Proteomics_normalized <- read_excel(
  xlsx_file,
  sheet = "Abundance_normalized"
)

########################################################################
################ Exploring the normalized proteomics data #############
########################################################################

# View the first few rows of the data
head(Proteomics_normalized)

# View the data
View(Proteomics_normalized)

# Check dimensions
dim(Proteomics_normalized)

# Check whether missing values are present
any(is.na(Proteomics_normalized))

# Count the total number of missing values
sum(is.na(Proteomics_normalized))


########################################################################
################ Formatting the proteomics data ########################
########################################################################

# Convert the proteomics data to a data frame
Proteomics_normalized <- as.data.frame(Proteomics_normalized)

# Use Protein_Accession as row names
rownames(Proteomics_normalized) <- Proteomics_normalized$Protein_Accession

# Remove the Protein_Accession column after assigning it as row names
Proteomics_normalized_2 <- Proteomics_normalized[, -1]

# Convert proteomics abundance values to numeric
Proteomics_normalized_2 <- data.frame(
  lapply(
    Proteomics_normalized_2,
    function(x) as.numeric(as.character(x))
  ),
  check.names = FALSE,
  row.names = rownames(Proteomics_normalized_2)
)

# Check dimensions
dim(Proteomics_normalized_2)

# Transpose the data so that samples are rows and proteins are columns
Proteomics_normalized_2_t <- as.data.frame(t(Proteomics_normalized_2))

# Check dimensions
dim(Proteomics_normalized_2_t)

# Check whether missing values are present before imputation
any(is.na(Proteomics_normalized_2_t))

# Count missing values before imputation
sum(is.na(Proteomics_normalized_2_t))

########################################################################
############ Checking completely missing proteins #####################
########################################################################

# Identify proteins for which all sample values are missing
all_na_proteins <- colSums(is.na(Proteomics_normalized_2_t)) ==
  nrow(Proteomics_normalized_2_t)

# Count proteins with missing values in all samples
sum(all_na_proteins)

# Check dimensions of the data used for imputation
dim(Proteomics_for_imputation)

########################################################################
################ missForest imputation ################################
########################################################################

# Set a random seed for reproducibility
# set.seed(123)

# Perform missing-value imputation using missForest
Proteomics_missForest <- missForest(
  Proteomics_normalized_2_t,
  maxiter = 10,
  ntree = 100,
  variablewise = FALSE,
  decreasing = FALSE,
  verbose = FALSE,
  mtry = floor(sqrt(ncol(Proteomics_normalized_2_t))),
  replace = TRUE,
  parallelize = "no"
)

# View the imputed data
View(Proteomics_missForest$ximp)

# Check the out-of-bag imputation error
Proteomics_missForest$OOBerror

########################################################################
################ Checking the imputation ###############################
########################################################################

# Extract the imputed data matrix
Proteomics_imputed <- Proteomics_missForest$ximp

# Count missing values before imputation
total_na_before <- sum(is.na(Proteomics_for_imputation))

# Count missing values after imputation
total_na_after <- sum(is.na(Proteomics_imputed))

# Display the number of missing values before and after imputation
cat("Missing values before imputation:", total_na_before, "\n")
cat("Missing values after imputation:", total_na_after, "\n")

########################################################################
################ Saving the imputed proteomics data ###################
########################################################################

# Add sample names as a column
imputed_data_with_samples <- cbind(
  Sample = rownames(Proteomics_imputed),
  Proteomics_imputed
)

# Save the imputed proteomics data as an Excel file
write_xlsx(
  imputed_data_with_samples,
  "Proteomics_imputed_data.xlsx"
)

########################################################################
######################## Saving R workspace ############################
########################################################################

# Save all objects in the current R session's global environment (workspace).

save.image(file = "Proteomics_Missing_Value_Imputation.RData")