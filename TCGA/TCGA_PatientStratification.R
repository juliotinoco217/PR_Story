library(readr)

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

# Filter for ER-positive, PR-positive, and HER2-negative patients only
filtered_data <- selected_data[selected_data$ER_STATUS_BY_IHC == "Positive" & 
                             selected_data$PR_STATUS_BY_IHC == "Positive" &
                             selected_data$IHC_HER2 == "Negative" &
                             selected_data$ER_STATUS_IHC_PERCENT_POSITIVE != "[Not Available]" &
                             selected_data$PR_STATUS_IHC_PERCENT_POSITIVE != "[Not Available]", ]

# Add PR_Quartile based on PR_STATUS_IHC_PERCENT_POSITIVE
filtered_data$PR_Quartile <- "Unknown"
filtered_data$PR_Quartile[filtered_data$PR_STATUS_IHC_PERCENT_POSITIVE %in% c("<10%", "10-19%", "20-29%")] <- "PR-Low"
filtered_data$PR_Quartile[filtered_data$PR_STATUS_IHC_PERCENT_POSITIVE %in% c("30-39%", "40-49%", "50-59%")] <- "PR-Medium"
filtered_data$PR_Quartile[filtered_data$PR_STATUS_IHC_PERCENT_POSITIVE %in% c("60-69%", "70-79%", "80-89%", "90-99%")] <- "PR-High"

# Write the filtered data to a new file
write.table(filtered_data, file="selected_clinical_data.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Print summary of the data
cat("Number of ER+/PR+/HER2- patients with available receptor status:", nrow(filtered_data), "\n")
cat("Number of columns:", ncol(filtered_data), "\n")

# Print summary of receptor status to verify filtering
cat("\nVerifying receptor status distribution:\n")
cat("ER status:\n")
print(table(filtered_data$ER_STATUS_BY_IHC))
cat("\nPR status:\n")
print(table(filtered_data$PR_STATUS_BY_IHC))
cat("\nHER2 status:\n")
print(table(filtered_data$IHC_HER2))

# Print PR quartile distribution
cat("\nPR Quartile distribution:\n")
print(table(filtered_data$PR_Quartile))

# Print detailed breakdown of PR percentages within each quartile
cat("\nDetailed PR percentage breakdown by quartile:\n")
print(table(filtered_data$PR_STATUS_IHC_PERCENT_POSITIVE, filtered_data$PR_Quartile))
