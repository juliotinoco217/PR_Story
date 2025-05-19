# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)

# Set working directory and read data
# Read PR stratification data
pr_data <- read.delim("../PR_stratificationIDs.txt")

# Read CIBERSORT results
cibersort_low <- read.delim("../CibersortResults_PR_Low .txt")
cibersort_mid <- read.delim("../CibersortResults_PRmid.txt")
cibersort_high <- read.delim("../CIBERSORTx_Results_PR_high.txt")

# Add PR category to each dataset
cibersort_low$PR_Category <- "Low"
cibersort_mid$PR_Category <- "Medium"
cibersort_high$PR_Category <- "High"

# Combine all datasets
cibersort_all <- bind_rows(cibersort_low, cibersort_mid, cibersort_high)

# Set PR_Category as factor with correct order
cibersort_all$PR_Category <- factor(cibersort_all$PR_Category, 
                                   levels = c("Low", "Medium", "High"))

# Define cell types
cell_types <- c(
  # B cells
  "B.cells.naive", "B.cells.memory", "Plasma.cells",
  
  # T cells
  "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting",
  "T.cells.CD4.memory.activated", "T.cells.follicular.helper",
  "T.cells.regulatory..Tregs.", "T.cells.gamma.delta",
  
  # NK cells
  "NK.cells.resting", "NK.cells.activated",
  
  # Myeloid cells
  "Monocytes", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2",
  "Dendritic.cells.resting", "Dendritic.cells.activated",
  
  # Other immune cells
  "Mast.cells.resting", "Mast.cells.activated",
  "Eosinophils", "Neutrophils"
)

# Create output directory for results
results_dir <- "METABRIC Statistics"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Function to perform Kruskal-Wallis and Dunn's test
perform_kw_dunn <- function(data, cell_types) {
  # Initialize dataframes for results
  kw_results <- data.frame(
    Cell_Type = character(),
    KW_Statistic = numeric(),
    KW_P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  dunn_results <- data.frame(
    Cell_Type = character(),
    Comparison = character(),
    Statistic = numeric(),
    P_Value = numeric(),
    P_Adjusted = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Perform tests for each cell type
  for (cell_type in cell_types) {
    # Prepare data for current cell type
    current_data <- data.frame(
      Value = data[[cell_type]],
      PR_Category = data$PR_Category
    )
    
    # Kruskal-Wallis test
    kw_test <- kruskal.test(Value ~ PR_Category, data = current_data)
    
    # Add KW results to summary
    kw_results <- rbind(
      kw_results,
      data.frame(
        Cell_Type = cell_type,
        KW_Statistic = kw_test$statistic,
        KW_P_Value = kw_test$p.value
      )
    )
    
    # Perform Dunn's test for all cell types
    dunn <- current_data %>%
      dunn_test(Value ~ PR_Category, p.adjust.method = "BH")
    
    # Add results to dunn_results
    dunn_results <- rbind(
      dunn_results,
      data.frame(
        Cell_Type = cell_type,
        Comparison = paste(dunn$group1, dunn$group2, sep="-"),
        Statistic = dunn$statistic,
        P_Value = dunn$p,
        P_Adjusted = dunn$p.adj
      )
    )
    
    # Create boxplot with significance
    p <- ggplot(current_data, aes(x = PR_Category, y = Value)) +
      geom_boxplot(aes(fill = PR_Category), alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
      scale_fill_manual(values = c("Low" = "#5C6B39", 
                                 "Medium" = "#2A7D8C",
                                 "High" = "#B22234")) +
      labs(title = paste(cell_type, "Distribution by PR Category"),
           subtitle = paste("Kruskal-Wallis p =", format.pval(kw_test$p.value, digits = 3)),
           y = "Cell Proportion",
           x = "PR Category") +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Save plot
    ggsave(
      filename = file.path(results_dir, paste0(cell_type, "_boxplot.pdf")),
      plot = p,
      width = 8,
      height = 6
    )
  }
  
  # Add FDR correction to KW results
  kw_results$KW_P_Adjusted <- p.adjust(kw_results$KW_P_Value, method = "BH")
  
  # Save results
  write.csv(kw_results, 
            file = file.path(results_dir, "kruskal_wallis_results.csv"),
            row.names = FALSE)
  
  write.csv(dunn_results,
            file = file.path(results_dir, "dunns_test_results.csv"),
            row.names = FALSE)
  
  # Create summary of significant differences
  significant_results <- dunn_results %>%
    filter(P_Adjusted < 0.05) %>%
    arrange(Cell_Type, P_Adjusted)
  
  write.csv(significant_results,
            file = file.path(results_dir, "significant_differences.csv"),
            row.names = FALSE)
  
  # Return all results
  return(list(
    kw_results = kw_results,
    dunn_results = dunn_results,
    significant_results = significant_results
  ))
}

# Perform Kruskal-Wallis and Dunn's tests
test_results <- perform_kw_dunn(cibersort_all, cell_types)

# Print summary of significant differences
cat("\nKruskal-Wallis test results:\n")
print(test_results$kw_results %>% 
        filter(KW_P_Adjusted < 0.05) %>% 
        arrange(KW_P_Adjusted))

cat("\nPairwise comparisons (Dunn's test) completed for all cell types.\n")
cat("Results saved in:", file.path(results_dir, "dunns_test_results.csv"), "\n")

if (nrow(test_results$significant_results) > 0) {
  cat("\nSignificant pairwise differences (Dunn's test):\n")
  print(test_results$significant_results)
  
  # Create summary of significant cell types
  significant_cell_types <- unique(test_results$significant_results$Cell_Type)
  cat("\nCell types with significant differences between PR categories:\n")
  print(significant_cell_types)
}

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), 
           file.path(results_dir, "session_info.txt")) 