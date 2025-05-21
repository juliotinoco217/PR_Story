# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining plots
library(readr)  # For read_csv

# Function to clean column names
clean_colnames <- function(names) {
  # Replace any variation of naive with a standard form
  names <- gsub("na.*ve_B_cells", "naive_B_cells", names, perl = TRUE)
  names <- gsub("naïve_B_cells", "naive_B_cells", names)
  names <- gsub("nave_B_cells", "naive_B_cells", names)
  return(names)
}

# Read and process Telepristone data
telepristone_xcell <- read_csv("deconvolution/xCell/Telepristone_GSE126765_xCell_results.csv",
                              show_col_types = FALSE,
                              locale = locale(encoding = "latin1")) %>%
  as.data.frame()
rownames(telepristone_xcell) <- telepristone_xcell$SampleID
telepristone_xcell <- telepristone_xcell[,-1]  # Remove SampleID column
colnames(telepristone_xcell) <- clean_colnames(colnames(telepristone_xcell))

# Add a column to identify the dataset source
telepristone_xcell$Dataset <- "Telepristone"

## Defining Cell Population Lists ##

# CD4 T cells
cd4_populations <- c(
  "CD4_naive_T_cells",
  "CD4_memory_T_cells",
  "CD4_T_cells",
  "CD4+ Tcm",
  "CD4+ Tem",
  "Th1_cells",
  "Th2_cells",
  "Tregs"
)

# CD8 T cells
cd8_populations <- c(
  "CD8_T_cells",
  "CD8_Tcm",
  "CD8_Tem",
  "CD8_naive_T_cells"
)

# NK cells
nk_populations <- c(
  "NK_cells",
  "NKT"
)

# B cells
b_populations <- c(
  "B_cells",
  "naive_B_cells",
  "pro_B_cells",
  "Memory_B_cells",
  "Class_switched_memory_B_cells",
  "Plasma_cells"
)

# Myeloid cells
myeloid_populations <- c(
  "Monocytes",
  "Macrophages",
  "Macrophages_M1",
  "Macrophages_M2",
  "DC",
  "cDC",
  "pDC",
  "iDC",
  "aDC"
)

# Granulocytes
granulocyte_populations <- c(
  "Neutrophils",
  "Basophils",
  "Eosinophils",
  "Mast_cells"
)

# Endothelial cells
endothelial_populations <- c(
  "Endothelial_cells",
  "ly_Endothelial_cells",
  "mv_Endothelial_cells"
)

# Epithelial cells
epithelial_populations <- c(
  "Epithelial_cells"
)

# Stromal cells
stromal_populations <- c(
  "Fibroblasts",
  "MSC",
  "Pericytes"
)

# Compiled Scores
compiled_scores <- c(
  "ImmuneScore",
  "StromaScore",
  "MicroenvironmentScore"
)

## Function to create cell type-specific dataframes ##
create_cell_df <- function(data, cell_populations, cell_type_name) {
  df <- data %>%
    select(all_of(c(cell_populations, "Dataset"))) %>%
    pivot_longer(cols = all_of(cell_populations),
                names_to = "Cell_Population",
                values_to = "Score") %>%
    mutate(Cell_Type = cell_type_name)
  return(df)
}

# Telepristone sample information
telepristone_samples <- list(
    GSM3612340 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612341 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612342 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "negative", menopause = "postmenopausal"),
    GSM3612343 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "negative", menopause = "postmenopausal"),
    GSM3612344 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "low positive", menopause = "premenopausal"),
    GSM3612345 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "low positive", menopause = "premenopausal"),
    GSM3612346 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612347 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612348 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612349 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612350 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612351 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612352 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "negative", menopause = "premenopausal"),
    GSM3612353 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "negative", menopause = "premenopausal"),
    GSM3612354 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612355 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612356 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612357 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612358 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "low positive", menopause = "postmenopausal"),
    GSM3612359 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "low positive", menopause = "postmenopausal"),
    GSM3612360 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612361 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612362 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612363 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612364 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612365 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612366 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "negative", pr_status = "negative", menopause = "postmenopausal"),
    GSM3612367 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "negative", pr_status = "negative", menopause = "postmenopausal"),
    GSM3612368 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612369 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612370 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612371 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612372 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "negative", pr_status = "negative", menopause = "premenopausal"),
    GSM3612373 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "negative", pr_status = "negative", menopause = "premenopausal"),
    GSM3612374 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612375 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612376 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612377 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612378 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612379 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612380 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612381 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612382 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612383 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612384 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612385 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612386 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612387 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612388 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612389 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612390 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "negative", menopause = "premenopausal"),
    GSM3612391 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "negative", menopause = "premenopausal"),
    GSM3612392 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612393 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612394 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612395 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612396 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612397 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612398 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612399 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612400 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612401 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612402 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612403 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612404 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612405 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612406 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612407 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612408 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612409 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612410 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612411 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612412 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612413 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612414 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612415 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612416 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612417 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612418 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612419 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612420 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612421 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612422 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612423 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612424 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612425 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612426 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612427 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612428 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612429 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612430 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612431 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612432 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612433 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612434 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612435 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612436 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612437 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612438 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612439 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612440 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612441 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612442 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612443 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612444 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612445 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612446 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612447 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612448 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612449 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612450 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612451 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612452 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612453 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612454 = list(specimen_type = "core needle biospy", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612455 = list(specimen_type = "surgical specimen", treatment = "TPA", er_status = "high positive", pr_status = "high positive", menopause = "premenopausal"),
    GSM3612456 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612457 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612458 = list(specimen_type = "core needle biospy", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal"),
    GSM3612459 = list(specimen_type = "surgical specimen", treatment = "placebo", er_status = "high positive", pr_status = "high positive", menopause = "postmenopausal")
)

# Create enhanced publication-ready theme
publication_theme <- theme_minimal() + theme(
  # Title styling
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
  # Axis styling
  axis.title = element_text(size = 12, face = "bold"),
  axis.text = element_text(size = 11, color = "black"),
  axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 10)),
  # Grid styling
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  # Border styling
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
  # Legend styling
  legend.title = element_text(size = 12, face = "bold"),
  legend.text = element_text(size = 11),
  # Margins and background
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  panel.background = element_rect(fill = "white"),
  strip.background = element_rect(fill = "white"),
  strip.text = element_text(size = 12, face = "bold")
)

# Function to create spaghetti plots
create_spaghetti_plot <- function(data, x_var, y_var, group_var, 
                                title = NULL, x_lab = NULL, y_lab = NULL,
                                color_palette = NULL) {
  
  # Base plot
  p <- ggplot(data, aes_string(x = x_var, y = y_var, group = group_var)) +
    # Add individual patient lines
    geom_line(color = "gray70", linewidth = 0.5, alpha = 0.6) +
    # Add points for pre and post with slightly larger size
    geom_point(aes_string(color = x_var), size = 3.5) +
    # Facet by cell type
    facet_wrap(~Cell_Type, scales = "free_y", ncol = 2) +
    # Ensure y-axis starts at 0 and includes all data points
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1)),  # No expansion below 0, 10% above max
      limits = function(x) c(0, max(x) * 1.1)  # Start at 0, extend 10% above max
    )
  
  # Add custom color palette if provided
  if (!is.null(color_palette)) {
    p <- p + scale_color_manual(values = color_palette)
  }
  
  # Add labels
  p <- p + labs(title = title,
                x = NULL,  # Remove x-axis label since it's clear from the legend
                y = y_lab,
                color = "Time Point")  # Legend title
  
  # Apply publication theme with additional styling
  p <- p + publication_theme + 
    theme(
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      legend.position = "top",
      # Adjust facet appearance
      strip.text = element_text(size = 11, face = "bold"),
      panel.spacing = unit(1, "lines")
    )
  
  return(p)
}

# Create plots directory if it doesn't exist
dir.create("deconvolution/xCell/plots", showWarnings = FALSE, recursive = TRUE)