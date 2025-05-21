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

# Read and process MIPRA data
mipra_xcell <- read_csv("deconvolution/xCell/MIPRA_xCell_results.csv", 
                       show_col_types = FALSE,
                       locale = locale(encoding = "latin1")) %>%
  as.data.frame()
rownames(mipra_xcell) <- mipra_xcell$SampleID
mipra_xcell <- mipra_xcell[,-1]  # Remove SampleID column
colnames(mipra_xcell) <- clean_colnames(colnames(mipra_xcell))

# Print column names to check available cell types
print("Available cell types in MIPRA dataset:")
print(colnames(mipra_xcell))

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

#Compiled Scores ## 
compiled_scores <- c(
  "ImmuneScore",
  "StromaScore",
  "MicroenvironmentScore"
)

## MIPRA Patient Information ##

# Create a list of patient IDs with their pre and post treatment samples
# and their response status
mipra_patient_samples <- list(
    M062 = list(samples = c(pre = "M062B", post = "M062S"), response = "Response"),
    M070 = list(samples = c(pre = "M070B", post = "M070S"), response = "Response"),
    M073 = list(samples = c(pre = "M073B", post = "M073S"), response = "Non-response"),
    M077 = list(samples = c(pre = "M077B", post = "M077S"), response = "Response"),
    M090 = list(samples = c(pre = "M090B", post = "M090S"), response = "Non-response"),
    M094 = list(samples = c(pre = "M094B", post = "M094S"), response = "Non-response"),
    M105 = list(samples = c(pre = "M105B", post = "M105S"), response = "Non-response"),
    M140 = list(samples = c(pre = "M140B", post = "M140S"), response = "Response")
)

#### MIPRA CD8 Spaghetti Plot ### 

# Define colors for pre and post (keeping Pre first in legend)
mipra_colors <- c(
  "Pre" = "#2A7D8C",   # Deep Teal
  "Post" = "#B22234"   # Crimson Red
)

# Initialize empty data frame for CD8 data
mipra_cd8_data <- data.frame()

# Debug: Print raw CD8_Tcm values before processing
print("Raw CD8_Tcm values:")
print(mipra_xcell[, "CD8_Tcm", drop=FALSE])

# Process CD8 data for all patients
for(patient in names(mipra_patient_samples)) {
  pre_sample <- mipra_patient_samples[[patient]]$samples["pre"]
  post_sample <- mipra_patient_samples[[patient]]$samples["post"]
  
  # Get CD8 data for pre samples
  pre_data <- mipra_xcell[pre_sample, cd8_populations] %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), 
                names_to = "Cell_Type", 
                values_to = "Proportion") %>%
    mutate(
      TimePoint = factor("Pre", levels = c("Pre", "Post")),
      Patient = patient,
      # Clean up the cell type names for display
      Cell_Type = case_when(
        Cell_Type == "CD8_T_cells" ~ "CD8+ T cells",
        Cell_Type == "CD8_naive_T_cells" ~ "CD8+ naive T cells",
        Cell_Type == "CD8_Tcm" ~ "CD8+ central memory T cells",
        Cell_Type == "CD8_Tem" ~ "CD8+ effector memory T cells",
        TRUE ~ Cell_Type
      )
    )
  
  # Get CD8 data for post samples
  post_data <- mipra_xcell[post_sample, cd8_populations] %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), 
                names_to = "Cell_Type", 
                values_to = "Proportion") %>%
    mutate(
      TimePoint = factor("Post", levels = c("Pre", "Post")),
      Patient = patient,
      # Clean up the cell type names for display
      Cell_Type = case_when(
        Cell_Type == "CD8_T_cells" ~ "CD8+ T cells",
        Cell_Type == "CD8_naive_T_cells" ~ "CD8+ naive T cells",
        Cell_Type == "CD8_Tcm" ~ "CD8+ central memory T cells",
        Cell_Type == "CD8_Tem" ~ "CD8+ effector memory T cells",
        TRUE ~ Cell_Type
      )
    )
  
  mipra_cd8_data <- rbind(mipra_cd8_data, pre_data, post_data)
}

# Debug: Print all CD8+ central memory T cells data
print("\nAll CD8+ central memory T cells data after processing:")
print(mipra_cd8_data %>% 
  filter(Cell_Type == "CD8+ central memory T cells") %>%
  arrange(TimePoint, Patient) %>%
  arrange(desc(Proportion)) %>%
  select(Patient, TimePoint, Proportion))

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

# Function to create spaghetti plot
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

# Create the spaghetti plot using the custom function
mipra_cd8_plot <- create_spaghetti_plot(
  data = mipra_cd8_data,
  x_var = "TimePoint",
  y_var = "Proportion",
  group_var = "Patient",
  title = "CD8+ T Cell Populations in MIPRA Study",
  y_lab = "Cell Proportion",
  color_palette = mipra_colors
)

# Create plots directory if it doesn't exist
dir.create("deconvolution/xCell/plots", showWarnings = FALSE, recursive = TRUE)

# Save the plot
ggsave(
  filename = "deconvolution/xCell/plots/MIPRA_CD8_spaghetti.pdf",
  plot = mipra_cd8_plot,
  width = 10,
  height = 8,  # Slightly taller to accommodate facets
  units = "in",
  dpi = 300)