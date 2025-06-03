########################################################
# 1) Load packages
########################################################
install.packages("survival")
install.packages("survminer")
install.packages("dplyr")
install.packages("skimr")  # for data quality checks

library(survival)
library(survminer)
library(dplyr)
library(skimr)

########################################################
# 2) Read clinical data and check data quality
########################################################
# Read clinical data
clinical_data <- read.table(
  "selected_clinical_data.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  fill = TRUE
)

# Print initial data quality check
print("Initial data quality check:")
print("Number of records before conversion:")
print(nrow(clinical_data))

# Check for any non-numeric values in month columns
print("\nChecking for non-numeric values in month columns:")
print("OS_MONTHS unique values:")
print(unique(clinical_data$OS_MONTHS))
print("\nDFS_MONTHS unique values:")
print(unique(clinical_data$DFS_MONTHS))

# Convert time variables to numeric and create event indicators
clinical_data <- clinical_data %>%
  mutate(across(ends_with("_MONTHS"), ~ as.numeric(.x)))

# Print summary of NA values after conversion
print("\nNA counts after numeric conversion:")
print(summarise(clinical_data, 
  across(ends_with("_MONTHS"), ~ sum(is.na(.)))
))

# Full data quality report
print("\nFull data quality report:")
print(skim(clinical_data))

# Create event indicators and censor data
clinical_data <- clinical_data %>%
  mutate(
    OS_event = ifelse(OS_STATUS == "1:DECEASED", 1, 0),
    RFS_event = ifelse(DFS_STATUS == "1:Recurred/Progressed", 1, 0)
  ) %>%
  # Properly censor data at 150 months
  mutate(
    OS_event = ifelse(OS_MONTHS > 150, 0, OS_event),
    OS_MONTHS = pmin(OS_MONTHS, 150),
    RFS_event = ifelse(DFS_MONTHS > 150, 0, RFS_event),
    DFS_MONTHS = pmin(DFS_MONTHS, 150)
  )

# Print final data quality check
print("\nFinal data quality check after censoring:")
print("Number of records:")
print(nrow(clinical_data))
print("\nSummary of time variables:")
print(summary(clinical_data[c("OS_MONTHS", "DFS_MONTHS")]))
print("\nSummary of event indicators:")
print(table(clinical_data$OS_event, useNA = "ifany"))
print(table(clinical_data$RFS_event, useNA = "ifany"))

########################################################
# 3) Three-way PR Group Analysis
########################################################
# Create PR grouping
merged_data_groups <- clinical_data %>%
  mutate(
    PR_Group = case_when(
      PR_Quartile == "PR-High" ~ "PR_High",
      PR_Quartile == "PR-Medium" ~ "PR_Medium",
      PR_Quartile == "PR-Low" ~ "PR_Low",
      TRUE ~ NA_character_
    )
  )

# Print summary of PR groups
print("PR Group Summary:")
print(table(merged_data_groups$PR_Group, useNA = "always"))

# Convert PR_Group to factor with specific order
merged_data_groups$PR_Group <- factor(
  merged_data_groups$PR_Group,
  levels = c("PR_High", "PR_Medium", "PR_Low")
)

# Remove NA values
merged_data_groups <- merged_data_groups %>%
  filter(!is.na(PR_Group))

# Overall Survival Analysis for Groups
surv_obj_groups <- Surv(
  time = merged_data_groups$OS_MONTHS,
  event = merged_data_groups$OS_event
)

fit_km_groups <- survfit(
  surv_obj_groups ~ PR_Group,
  data = merged_data_groups
)

# Create OS plot for groups
p_groups_os <- ggsurvplot(
  fit_km_groups,
  data = merged_data_groups,
  conf.int = FALSE,
  risk.table = TRUE,
  pval = TRUE,
  palette = c("#2A7D8C", "#5C6B39", "#B22234"),  # High = blue, Medium = green, Low = red
  title = "Overall Survival by PR Expression Groups",
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  risk.table.height = 0.3,
  legend.labs = c("PR High", "PR Medium", "PR Low"),
  legend.title = "PR Status",
  xlim = c(0, 150)
)

# Save OS plot
png("./TCGA_Plot_OS.png", width = 10, height = 8, units = "in", res = 300)
print(p_groups_os)
dev.off()

pdf("./TCGA_Plot_OS.pdf", width = 10, height = 8)
print(p_groups_os)
dev.off()

# RFS Analysis for Groups
surv_obj_rfs_groups <- Surv(
  time = merged_data_groups$DFS_MONTHS,
  event = merged_data_groups$RFS_event
)

fit_km_rfs_groups <- survfit(
  surv_obj_rfs_groups ~ PR_Group,
  data = merged_data_groups
)

# Create RFS plot for groups
p_groups_rfs <- ggsurvplot(
  fit_km_rfs_groups,
  data = merged_data_groups,
  conf.int = FALSE,
  risk.table = TRUE,
  pval = TRUE,
  palette = c("#2A7D8C", "#5C6B39", "#B22234"),  # High = blue, Medium = green, Low = red
  title = "Relapse Free Survival by PR Expression Groups",
  xlab = "Time (months)",
  ylab = "Relapse Free Survival Probability",
  risk.table.height = 0.3,
  legend.labs = c("PR High", "PR Medium", "PR Low"),
  legend.title = "PR Status",
  xlim = c(0, 150)
)

# Save RFS plot
png("./TCGA_Plot_RFS.png", width = 10, height = 8, units = "in", res = 300)
print(p_groups_rfs)
dev.off()

pdf("./TCGA_Plot_RFS.pdf", width = 10, height = 8)
print(p_groups_rfs)
dev.off()

# Print summary statistics for group analyses
# Open file connection to save statistics
sink("TCGA_survival_statistics.txt")

print("=================================================================")
print("SURVIVAL ANALYSIS RESULTS - PR Expression Groups")
print("=================================================================")

print("\nGroup Sizes:")
print(table(merged_data_groups$PR_Group))

print("\n=================================================================")
print("THREE-WAY COMPARISON")
print("=================================================================")

print("\nOverall Survival Analysis:")
print("---------------------------")
print(summary(fit_km_groups))
print("\nLog-rank test results:")
print(survdiff(surv_obj_groups ~ PR_Group, data = merged_data_groups))

print("\nRelapse Free Survival Analysis:")
print("--------------------------------")
print(summary(fit_km_rfs_groups))
print("\nLog-rank test results:")
print(survdiff(surv_obj_rfs_groups ~ PR_Group, data = merged_data_groups))

print("\n=================================================================")
print("PAIRWISE COMPARISONS")
print("=================================================================")

print("\nPairwise Comparisons - Overall Survival:")
print("----------------------------------------")
print("\nHigh vs Medium:")
high_med <- merged_data_groups %>% filter(PR_Group %in% c("PR_High", "PR_Medium"))
print(survdiff(Surv(OS_MONTHS, OS_event) ~ PR_Group, data = high_med))

print("\nHigh vs Low:")
high_low <- merged_data_groups %>% filter(PR_Group %in% c("PR_High", "PR_Low"))
print(survdiff(Surv(OS_MONTHS, OS_event) ~ PR_Group, data = high_low))

print("\nMedium vs Low:")
med_low <- merged_data_groups %>% filter(PR_Group %in% c("PR_Medium", "PR_Low"))
print(survdiff(Surv(OS_MONTHS, OS_event) ~ PR_Group, data = med_low))

print("\nPairwise Comparisons - Relapse Free Survival:")
print("---------------------------------------------")
print("\nHigh vs Medium:")
print(survdiff(Surv(DFS_MONTHS, RFS_event) ~ PR_Group, data = high_med))

print("\nHigh vs Low:")
print(survdiff(Surv(DFS_MONTHS, RFS_event) ~ PR_Group, data = high_low))

print("\nMedium vs Low:")
print(survdiff(Surv(DFS_MONTHS, RFS_event) ~ PR_Group, data = med_low))

# Close the file connection 