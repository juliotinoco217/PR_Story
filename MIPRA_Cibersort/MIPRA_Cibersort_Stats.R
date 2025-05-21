# Load required libraries
required_packages <- c("ggplot2", "dplyr", "tidyr")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load packages
for(pkg in required_packages) {
  suppressMessages(library(pkg, character.only = TRUE))
}

# Create output directory for results
results_dir <- "MIPRA_Statistics"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Read CIBERSORT results
cibersort_results <- read.delim("deconvolution/CibersortxResults_MIPRA_.txt", row.names = 1) %>%
  mutate(
    Sample_ID = rownames(.),
    Treatment = factor(ifelse(grepl("B$", Sample_ID), "Pre", "Post"), levels = c("Pre", "Post")),
    Patient = sub("[BS]$", "", Sample_ID)
  ) %>%
  group_by(Patient) %>%
  filter(n() == 2) %>%  # Keep only complete pairs
  ungroup()

# Define patient response information
patient_samples <- list(
    M062 = list(samples = c(pre = "M062B", post = "M062S"), response = "Response"),
    M070 = list(samples = c(pre = "M070B", post = "M070S"), response = "Response"),
    M073 = list(samples = c(pre = "M073B", post = "M073S"), response = "Non-response"),
    M077 = list(samples = c(pre = "M077B", post = "M077S"), response = "Response"),
    M090 = list(samples = c(pre = "M090B", post = "M090S"), response = "Non-response"),
    M094 = list(samples = c(pre = "M094B", post = "M094S"), response = "Non-response"),
    M105 = list(samples = c(pre = "M105B", post = "M105S"), response = "Non-response"),
    M140 = list(samples = c(pre = "M140B", post = "M140S"), response = "Response")
)

# Create response mapping
response_mapping <- data.frame(
  Patient = names(patient_samples),
  Response = sapply(patient_samples, function(x) x$response)
)

# Add response information
cibersort_with_response <- cibersort_results %>%
  left_join(response_mapping, by = "Patient")

# Define cell types
cell_types <- c(
  "B cells naive", "B cells memory", "Plasma cells",
  "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta",
  "NK cells resting", "NK cells activated",
  "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
  "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated",
  "Eosinophils", "Neutrophils"
)

# Function to perform paired t-tests
perform_paired_tests <- function(data, cell_types) {
  test_results <- data.frame(
    Cell_Type = character(),
    Statistic = numeric(),
    P_Value = numeric(),
    Mean_Diff = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cell_type in cell_types) {
    # Prepare data for current cell type
    test_data <- data %>%
      select(Patient, Treatment, !!sym(cell_type)) %>%
      spread(Treatment, !!sym(cell_type))
    
    # Calculate differences
    differences <- test_data$Post - test_data$Pre
    
    # Perform paired t-test
    t_test <- tryCatch({
      test <- t.test(test_data$Post, test_data$Pre, paired = TRUE)
      
      data.frame(
        Cell_Type = cell_type,
        Statistic = test$statistic,
        P_Value = test$p.value,
        Mean_Diff = mean(differences),
        CI_Lower = test$conf.int[1],
        CI_Upper = test$conf.int[2]
      )
    }, error = function(e) NULL)
    
    if (!is.null(t_test)) {
      test_results <- rbind(test_results, t_test)
    }
  }
  
  # Add FDR correction
  test_results$P_Adjusted <- p.adjust(test_results$P_Value, method = "BH")
  return(test_results)
}

# Function to perform response group analysis
perform_response_analysis <- function(data, cell_types) {
  response_results <- data.frame(
    Cell_Type = character(),
    Timepoint = character(),
    Statistic = numeric(),
    P_Value = numeric(),
    Mean_Diff = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cell_type in cell_types) {
    for (time_point in c("Pre", "Post")) {
      current_data <- data %>%
        filter(Treatment == time_point)
      
      # Perform t-test
      test_result <- tryCatch({
        response_data <- split(current_data[[cell_type]], current_data$Response)
        test <- t.test(response_data[["Response"]], response_data[["Non-response"]])
        
        data.frame(
          Cell_Type = cell_type,
          Timepoint = time_point,
          Statistic = test$statistic,
          P_Value = test$p.value,
          Mean_Diff = diff(test$estimate),
          CI_Lower = test$conf.int[1],
          CI_Upper = test$conf.int[2]
        )
      }, error = function(e) NULL)
      
      if (!is.null(test_result)) {
        response_results <- rbind(response_results, test_result)
      }
    }
  }
  
  # Add FDR correction
  response_results$P_Adjusted <- p.adjust(response_results$P_Value, method = "BH")
  return(response_results)
}

# Perform analyses
paired_results <- perform_paired_tests(cibersort_results, cell_types)
response_results <- perform_response_analysis(cibersort_with_response, cell_types)

# Save results
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
write.csv(paired_results, file.path(results_dir, paste0("paired_ttest_results_", timestamp, ".csv")), row.names = FALSE)
write.csv(response_results, file.path(results_dir, paste0("response_ttest_results_", timestamp, ".csv")), row.names = FALSE)

# Print significant results (P < 0.05)
significant_paired <- paired_results %>% filter(P_Value < 0.05)
significant_response <- response_results %>% filter(P_Value < 0.05)

if(nrow(significant_paired) > 0) print(significant_paired)
if(nrow(significant_response) > 0) print(significant_response)
