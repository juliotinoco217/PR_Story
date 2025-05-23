# Read the clinical data file
data <- read.delim("brca_tcga/data_clinical_patient.txt", header=TRUE, sep="\t", comment.char="#")

# Select the relevant columns
selected_columns <- c(
  "PATIENT_ID",
  "ER_STATUS_BY_IHC",
  "ER_STATUS_IHC_PERCENT_POSITIVE",
  "PR_STATUS_BY_IHC", 
  "PR_STATUS_IHC_PERCENT_POSITIVE",
  "IHC_HER2",
  "HER2_FISH_STATUS",
  "OS_STATUS",
  "OS_MONTHS",
  "DFS_STATUS",
  "DFS_MONTHS"
)

# Create a new dataframe with selected columns
selected_data <- data[, selected_columns]

# Write the selected data to a new file
write.table(selected_data, file="selected_clinical_data.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Print summary of the data
cat("Number of patients:", nrow(selected_data), "\n")
cat("Number of columns:", ncol(selected_data), "\n")
