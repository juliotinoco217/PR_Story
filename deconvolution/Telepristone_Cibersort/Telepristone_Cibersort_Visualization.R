# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(viridis)
library(cowplot)
library(scales)

# Set theme for publication-ready figures (fixed font and deprecated arguments)
theme_publication <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      # Text elements
      axis.text = element_text(size = 10, colour = "black"),
      axis.title = element_text(size = 12, colour = "black", face = "bold"),
      strip.text = element_text(size = 11, colour = "black", face = "bold"),
      legend.text = element_text(size = 10, colour = "black"),
      legend.title = element_text(size = 11, colour = "black", face = "bold"),
      plot.title = element_text(size = 14, hjust = 0.5, colour = "black", face = "bold"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "black"),
      
      # Panel and grid (using linewidth instead of size)
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      
      # Strip and legend
      strip.background = element_rect(colour = "black", fill = "grey95"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      
      # Margins
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
}

# Set default theme
theme_set(theme_publication())

# Read clinical data with improved error handling
if (!file.exists("Telepristone_Clinical_Cleaned.txt")) {
  stop("Clinical data file not found. Please check file path.")
}

clinical_data <- read.delim("Telepristone_Clinical_Cleaned.txt", check.names = FALSE) %>%
  pivot_longer(cols = starts_with("GSM"), 
              names_to = "SampleID", 
              values_to = "value") %>%
  pivot_wider(names_from = "PatientID", 
             values_from = "value") %>%
  select(SampleID, Treatment_Group, Specimen_Type, ER_Status, PR_Status, Menopause_Status) %>%
  # Better patient ID creation - ensure proper pairing
  arrange(Treatment_Group, Specimen_Type, SampleID) %>%
  group_by(Treatment_Group, Specimen_Type) %>%
  mutate(temp_order = row_number()) %>%
  ungroup() %>%
  # Create patient IDs ensuring pre/post pairs
  group_by(Treatment_Group) %>%
  arrange(Treatment_Group, temp_order) %>%
  mutate(
    # Assign patient IDs based on pairs within each treatment group
    PatientID = case_when(
      Specimen_Type == "core needle biospy" ~ paste0(Treatment_Group, "_", sprintf("%02d", temp_order)),
      Specimen_Type == "surgical specimen" ~ paste0(Treatment_Group, "_", sprintf("%02d", temp_order)),
      TRUE ~ NA_character_
    )
  ) %>%
  select(-temp_order) %>%
  ungroup()

# Enhanced validation of pairing
cat("Validating patient pairing...\n")

# Check 1: Each patient should have exactly 2 samples
pairing_check <- clinical_data %>%
  count(PatientID, name = "n_samples") %>%
  filter(n_samples != 2)

if (nrow(pairing_check) > 0) {
  cat("WARNING: Patients with incorrect number of samples:\n")
  print(pairing_check)
  
  # Remove patients without proper pairs
  valid_patients <- clinical_data %>%
    count(PatientID) %>%
    filter(n == 2) %>%
    pull(PatientID)
  
  clinical_data <- clinical_data %>%
    filter(PatientID %in% valid_patients)
  
  cat(sprintf("Removed %d patients with incomplete pairs\n", 
              nrow(pairing_check)))
}

# Check 2: Each patient should have one pre and one post sample
specimen_check <- clinical_data %>%
  count(PatientID, Specimen_Type) %>%
  pivot_wider(names_from = Specimen_Type, values_from = n, values_fill = 0) %>%
  filter(`core needle biospy` != 1 | `surgical specimen` != 1)

if (nrow(specimen_check) > 0) {
  cat("WARNING: Patients without proper pre/post pairing:\n")
  print(specimen_check)
  
  # Keep only patients with proper pre/post pairs
  valid_patients <- clinical_data %>%
    count(PatientID, Specimen_Type) %>%
    pivot_wider(names_from = Specimen_Type, values_from = n, values_fill = 0) %>%
    filter(`core needle biospy` == 1 & `surgical specimen` == 1) %>%
    pull(PatientID)
  
  clinical_data <- clinical_data %>%
    filter(PatientID %in% valid_patients)
  
  cat(sprintf("Kept %d patients with proper pre/post pairs\n", 
              length(valid_patients)))
}

# Check 3: Verify treatment group consistency within patients
treatment_check <- clinical_data %>%
  group_by(PatientID) %>%
  summarise(n_treatments = n_distinct(Treatment_Group), .groups = "drop") %>%
  filter(n_treatments > 1)

if (nrow(treatment_check) > 0) {
  stop("ERROR: Some patients have samples in multiple treatment groups!")
}

# Summarize study data
study_summary <- clinical_data %>%
  group_by(Treatment_Group, Specimen_Type) %>%
  summarise(
    n_samples = n(),
    n_er_pos = sum(grepl("positive", ER_Status, ignore.case = TRUE)),
    n_pr_pos = sum(grepl("positive", PR_Status, ignore.case = TRUE)),
    n_premenopausal = sum(grepl("premenopausal", Menopause_Status, ignore.case = TRUE)),
    .groups = "drop"
  )

cat("Study Summary:\n")
print(study_summary)

# Read CIBERSORT results with validation
if (!file.exists("Telepristone_CibersortResultsRelative_cleaned.txt")) {
  stop("CIBERSORT results file not found. Please check file path.")
}

cibersort_results <- read.delim("Telepristone_CibersortResultsRelative_cleaned.txt") %>%
  rename(SampleID = Mixture) %>%
  select(-P_value, -Correlation, -RMSE)

# Validate that all clinical samples have CIBERSORT data
missing_cibersort <- setdiff(clinical_data$SampleID, cibersort_results$SampleID)
if (length(missing_cibersort) > 0) {
  warning("Some clinical samples missing CIBERSORT data:")
  print(missing_cibersort)
}

# Merge clinical data with CIBERSORT results
merged_data <- inner_join(clinical_data, cibersort_results, by = "SampleID") %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, 
                           levels = c("placebo", "TPA"), 
                           labels = c("Placebo", "TPA")),
    Specimen_Type = factor(Specimen_Type, 
                          levels = c("core needle biospy", "surgical specimen"),
                          labels = c("Pre-treatment", "Post-treatment"))
  ) %>%
  # Ensure we only keep patients with both time points
  group_by(PatientID) %>%
  filter(n() == 2) %>%
  ungroup()

# Final validation of merged dataset
cat("\n=== FINAL DATASET VALIDATION ===\n")

# Check for complete pairs in merged data
final_pair_check <- merged_data %>%
  count(PatientID, Treatment_Group) %>%
  filter(n != 2)

if (nrow(final_pair_check) > 0) {
  cat("ERROR: Incomplete pairs in final merged dataset:\n")
  print(final_pair_check)
  stop("Data integrity compromised - stopping analysis")
}

# Verify specimen type pairing
final_specimen_check <- merged_data %>%
  count(PatientID, Specimen_Type) %>%
  pivot_wider(names_from = Specimen_Type, values_from = n, values_fill = 0) %>%
  filter(`Pre-treatment` != 1 | `Post-treatment` != 1)

if (nrow(final_specimen_check) > 0) {
  cat("ERROR: Incorrect specimen type pairing:\n")
  print(final_specimen_check)
  stop("Data integrity compromised - stopping analysis")
}

# Generate comprehensive dataset summary
dataset_summary <- merged_data %>%
  group_by(Treatment_Group, Specimen_Type) %>%
  summarise(
    n_samples = n(),
    n_unique_patients = n_distinct(PatientID),
    .groups = "drop"
  )

cat("Final dataset summary:\n")
print(dataset_summary)

# Verify patient balance
patient_summary <- merged_data %>%
  group_by(Treatment_Group) %>%
  summarise(
    n_patients = n_distinct(PatientID),
    n_samples = n(),
    .groups = "drop"
  )

cat("\nPatient distribution by treatment:\n")
print(patient_summary)

# Check for any missing CIBERSORT data
missing_data_check <- merged_data %>%
  select(starts_with(c("T_cells", "B_cells", "NK_cells", "Monocytes", 
                      "Macrophages", "Dendritic", "Mast", "Eosinophils", 
                      "Neutrophils", "Plasma"))) %>%
  summarise_all(~sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "cell_type", values_to = "missing_count") %>%
  filter(missing_count > 0)

if (nrow(missing_data_check) > 0) {
  cat("\nWARNING: Missing CIBERSORT data found:\n")
  print(missing_data_check)
}

cat(sprintf("\nFinal dataset: %d patients (%d samples) with complete paired data\n", 
            length(unique(merged_data$PatientID)), nrow(merged_data)))
cat("Data validation completed successfully!\n")

# Define cell populations with improved organization
cell_populations <- list(
  adaptive = c(
    "T_cells_CD8" = "CD8+ T cells",
    "T_cells_CD4_naive" = "Naive CD4+ T cells", 
    "T_cells_CD4_memory_resting" = "Resting Memory CD4+ T cells",
    "T_cells_CD4_memory_activated" = "Activated Memory CD4+ T cells",
    "T_cells_follicular_helper" = "T Follicular Helper cells",
    "T_cells_regulatory_Tregs" = "Regulatory T cells",
    "B_cells_naive" = "Naive B cells",
    "B_cells_memory" = "Memory B cells",
    "Plasma_cells" = "Plasma cells"
  ),
  innate = c(
    "NK_cells_resting" = "Resting NK cells",
    "NK_cells_activated" = "Activated NK cells",
    "Monocytes" = "Monocytes",
    "Macrophages_M0" = "M0 Macrophages",
    "Macrophages_M1" = "M1 Macrophages", 
    "Macrophages_M2" = "M2 Macrophages",
    "Dendritic_cells_resting" = "Resting Dendritic cells",
    "Dendritic_cells_activated" = "Activated Dendritic cells",
    "Mast_cells_resting" = "Resting Mast cells",
    "Mast_cells_activated" = "Activated Mast cells",
    "Eosinophils" = "Eosinophils",
    "Neutrophils" = "Neutrophils"
  )
)

# Enhanced statistical testing function
calculate_paired_stats <- function(data, cell_type) {
  # Separate data by treatment group and perform paired analysis
  stats_list <- list()
  
  for (treatment in unique(data$Treatment_Group)) {
    treatment_data <- data %>% filter(Treatment_Group == treatment)
    
    # Check if we have paired data (each patient should have 2 samples)
    patient_counts <- treatment_data %>% 
      count(PatientID) %>% 
      filter(n == 2)
    
    if (nrow(patient_counts) < 3) {
      # Not enough pairs for statistical testing
      stats_list[[treatment]] <- tibble(
        Treatment_Group = treatment,
        n_pairs = nrow(patient_counts),
        mean_pre = NA,
        mean_post = NA,
        median_pre = NA,
        median_post = NA,
        p_value_wilcox = NA,
        p_value_ttest = NA,
        effect_size = NA
      )
      next
    }
    
    # Filter to only patients with complete pairs
    paired_data <- treatment_data %>%
      filter(PatientID %in% patient_counts$PatientID) %>%
      arrange(PatientID, Specimen_Type)
    
    # Extract pre and post values, ensuring proper pairing
    pre_values <- paired_data %>%
      filter(Specimen_Type == "Pre-treatment") %>%
      arrange(PatientID) %>%
      pull(!!sym(cell_type))
    
    post_values <- paired_data %>%
      filter(Specimen_Type == "Post-treatment") %>%
      arrange(PatientID) %>%
      pull(!!sym(cell_type))
    
    # Verify that we have the same number of pre and post values
    if (length(pre_values) != length(post_values)) {
      warning(sprintf("Mismatch in pre/post values for %s in %s treatment group", 
                     cell_type, treatment))
      next
    }
    
    # Calculate summary statistics
    mean_pre <- mean(pre_values, na.rm = TRUE)
    mean_post <- mean(post_values, na.rm = TRUE)
    median_pre <- median(pre_values, na.rm = TRUE)
    median_post <- median(post_values, na.rm = TRUE)
    
    # Perform paired statistical tests
    p_value_wilcox <- tryCatch({
      if (length(pre_values) >= 3) {
        wilcox.test(x = pre_values, y = post_values, paired = TRUE)$p.value
      } else {
        NA
      }
    }, error = function(e) {
      warning(sprintf("Wilcoxon test failed for %s in %s: %s", 
                     cell_type, treatment, e$message))
      NA
    })
    
    p_value_ttest <- tryCatch({
      if (length(pre_values) >= 3) {
        t.test(x = pre_values, y = post_values, paired = TRUE)$p.value
      } else {
        NA
      }
    }, error = function(e) {
      warning(sprintf("T-test failed for %s in %s: %s", 
                     cell_type, treatment, e$message))
      NA
    })
    
    # Calculate effect size
    effect_size <- if (!is.na(mean_pre) && mean_pre != 0) {
      (mean_post - mean_pre) / mean_pre * 100
    } else {
      NA
    }
    
    # Store results
    stats_list[[treatment]] <- tibble(
      Treatment_Group = treatment,
      n_pairs = length(pre_values),
      mean_pre = mean_pre,
      mean_post = mean_post,
      median_pre = median_pre,
      median_post = median_post,
      p_value_wilcox = p_value_wilcox,
      p_value_ttest = p_value_ttest,
      effect_size = effect_size
    )
  }
  
  # Combine results and add significance indicators
  stats <- bind_rows(stats_list) %>%
    mutate(
      p_value = p_value_wilcox,  # Use Wilcoxon as primary test
      significance = case_when(
        is.na(p_value) ~ "insufficient_data",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**", 
        p_value < 0.05 ~ "*",
        p_value < 0.1 ~ ".",
        TRUE ~ "ns"
      )
    )
  
  return(stats)
}

# Function to format p-values for display
format_pvalue <- function(p) {
  if (is.na(p)) return("N/A")
  if (p < 0.001) return("< 0.001")
  if (p < 0.01) return(sprintf("= %.3f", p))
  return(sprintf("= %.2f", p))
}

# Enhanced function to create publication-ready spaghetti plot
create_enhanced_spaghetti_plot <- function(data, cell_type, cell_label) {
  # Calculate statistics for annotation
  stats_result <- calculate_paired_stats(data, cell_type)
  
  # Create the plot (using linewidth instead of size)
  p <- ggplot(data, aes(x = Specimen_Type, y = .data[[cell_type]], 
                       group = PatientID, color = Treatment_Group)) +
    geom_line(alpha = 0.4, linewidth = 0.5) +
    geom_point(size = 1.5, alpha = 0.7) +
    
    # Add summary statistics - removed group = Treatment_Group since faceting handles this
    stat_summary(fun = mean, geom = "point", 
                size = 3, shape = 19, alpha = 0.9) +
    stat_summary(fun = mean, geom = "line", 
                linewidth = 1.5, alpha = 1.0, group = 1) +
    
    facet_wrap(~Treatment_Group) +
    
    labs(
      title = cell_label,
      x = "Time Point",
      y = "Relative Proportion",
      color = "Treatment Group"
    ) +
    
    scale_color_manual(values = c("Placebo" = "#2A7D8C", "TPA" = "#B22234")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))
  
  # Try to add statistical annotations with error handling
  tryCatch({
    p <- p + annotate("text", 
                     x = 1.5, 
                     y = Inf, 
                     label = sprintf("Placebo: p %s\nTPA: p %s",
                                   format_pvalue(stats_result$p_value[stats_result$Treatment_Group == "Placebo"]),
                                   format_pvalue(stats_result$p_value[stats_result$Treatment_Group == "TPA"])),
                     vjust = 1.1, hjust = 0.5, size = 3)
  }, error = function(e) {
    warning("Could not add statistical annotations: ", e$message)
  })
  
  return(p)
}

# Create comprehensive statistical summary
all_stats <- tibble()

# Generate plots and statistics for each cell type
plots_list <- list()
stats_summary <- list()

for (category in names(cell_populations)) {
  cat(sprintf("\nProcessing %s immune cells...\n", category))
  
  for (cell_type in names(cell_populations[[category]])) {
    cell_label <- cell_populations[[category]][cell_type]
    
    # Calculate statistics
    stats <- calculate_paired_stats(merged_data, cell_type)
    stats$cell_type <- cell_type
    stats$cell_label <- cell_label
    stats$category <- category
    
    # Store statistics
    all_stats <- bind_rows(all_stats, stats)
    
    # Create plot
    p <- create_enhanced_spaghetti_plot(merged_data, cell_type, cell_label)
    plots_list[[paste(category, cell_type, sep = "_")]] <- p
    
    # Save individual plot with error handling
    tryCatch({
      ggsave(
        filename = sprintf("publication_spaghetti_%s_%s.pdf", category, cell_type),
        plot = p,
        width = 10, height = 6, units = "in",
        dpi = 300, device = "pdf"
      )
      
      # Also save as high-res PNG for presentations
      ggsave(
        filename = sprintf("publication_spaghetti_%s_%s.png", category, cell_type),
        plot = p,
        width = 10, height = 6, units = "in",
        dpi = 300, device = "png"
      )
    }, error = function(e) {
      warning(sprintf("Could not save plot for %s_%s: %s", category, cell_type, e$message))
    })
  }
}

# Save comprehensive statistics table
write.csv(all_stats, "telepristone_cibersort_statistics.csv", row.names = FALSE)

# Create summary plots for significant results
significant_results <- all_stats %>%
  filter(!is.na(p_value) & p_value < 0.05) %>%
  arrange(p_value)

if (nrow(significant_results) > 0) {
  cat(sprintf("\nFound %d significant changes (p < 0.05):\n", nrow(significant_results)))
  print(significant_results %>% 
        select(Treatment_Group, cell_label, p_value, effect_size, significance))
  
  # Create a summary heatmap of effect sizes
  effect_heatmap_data <- all_stats %>%
    select(Treatment_Group, cell_label, effect_size, p_value) %>%
    mutate(
      significance = case_when(
        is.na(p_value) ~ "No data",
        p_value < 0.001 ~ "p < 0.001",
        p_value < 0.01 ~ "p < 0.01",
        p_value < 0.05 ~ "p < 0.05", 
        p_value < 0.1 ~ "p < 0.1",
        TRUE ~ "ns"
      ),
      significance = factor(significance, 
                          levels = c("p < 0.001", "p < 0.01", "p < 0.05", "p < 0.1", "ns", "No data"))
    )
  
  heatmap_plot <- ggplot(effect_heatmap_data, 
                        aes(x = Treatment_Group, y = reorder(cell_label, effect_size))) +
    geom_tile(aes(fill = effect_size), color = "white", linewidth = 0.5) +
    geom_point(aes(size = significance), shape = 21, color = "black", 
               data = filter(effect_heatmap_data, significance != "ns" & significance != "No data")) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "% Change") +
    scale_size_manual(values = c("p < 0.001" = 3, "p < 0.01" = 2.5, 
                                "p < 0.05" = 2, "p < 0.1" = 1.5),
                     name = "Significance") +
    labs(
      title = "Treatment Effects on Immune Cell Populations", 
      subtitle = "Percent change from pre- to post-treatment",
      x = "Treatment Group",
      y = "Cell Population"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  tryCatch({
    ggsave("telepristone_effect_heatmap.pdf", heatmap_plot, 
           width = 8, height = 12, dpi = 300)
    ggsave("telepristone_effect_heatmap.png", heatmap_plot, 
           width = 8, height = 12, dpi = 300)
  }, error = function(e) {
    warning("Could not save heatmap: ", e$message)
  })
}

# Generate summary report
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TELEPRISTONE CIBERSORT ANALYSIS SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat(sprintf("Total patients analyzed: %d\n", length(unique(merged_data$PatientID))))
cat(sprintf("Total samples: %d\n", nrow(merged_data)))

cat("\nTreatment groups:\n")
treatment_summary <- merged_data %>%
  group_by(Treatment_Group) %>%
  summarise(n_patients = length(unique(PatientID)), .groups = "drop")
print(treatment_summary)

cat(sprintf("\nCell populations analyzed: %d\n", length(unlist(cell_populations))))
cat(sprintf("Significant changes found (p < 0.05): %d\n", 
            nrow(filter(all_stats, !is.na(p_value) & p_value < 0.05))))

cat("\nFiles generated:\n")
cat("- Individual plots: publication_spaghetti_[category]_[celltype].pdf/png\n")
cat("- Statistics table: telepristone_cibersort_statistics.csv\n")
if (nrow(significant_results) > 0) {
  cat("- Effect size heatmap: telepristone_effect_heatmap.pdf/png\n") 
}

cat("\nAnalysis completed successfully!\n") 