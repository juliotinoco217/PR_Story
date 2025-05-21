# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining plots

# Set the output directory for plots
plot_dir <- "MIPRA_Cibersortplots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Define the single_plot_theme
single_plot_theme <- theme_minimal() +
  theme(
    # Title styling
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    # Axis styling
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 10)),
    # Grid styling
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    # Border styling
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    # Background
    panel.background = element_rect(fill = "white"),
    # Margins
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

# Define patient samples
patient_samples <- list(
  M062 = list(pre = "M062B", post = "M062S"),
  M070 = list(pre = "M070B", post = "M070S"),
  M073 = list(pre = "M073B", post = "M073S"),
  M077 = list(pre = "M077B", post = "M077S"),
  M090 = list(pre = "M090B", post = "M090S"),
  M094 = list(pre = "M094B", post = "M094S"),
  M105 = list(pre = "M105B", post = "M105S"),
  M140 = list(pre = "M140B", post = "M140S")
)

# Read CIBERSORT results
cibersort_results <- read.delim("deconvolution/CibersortxResults_MIPRA_.txt", row.names = 1)

# Print initial data check
print("Initial data dimensions:")
print(dim(cibersort_results))
print("Column names:")
print(colnames(cibersort_results))

# Add sample IDs column from row names
cibersort_results$Sample_ID <- rownames(cibersort_results)

# Extract patient ID and treatment status
cibersort_results <- cibersort_results %>%
  mutate(
    Treatment = ifelse(grepl("B$", Sample_ID), "Pre", "Post"),
    Patient = sub("[BS]$", "", Sample_ID),
    Treatment = factor(Treatment, levels = c("Pre", "Post"))
  )

# Print data after processing
print("Data after processing:")
print(head(cibersort_results))

# Ensure we have complete pairs
complete_pairs <- cibersort_results %>%
  group_by(Patient) %>%
  filter(n() == 2) %>%  # Keep only patients with both pre and post samples
  ungroup()

# Print complete pairs info
print("Complete pairs dimensions:")
print(dim(complete_pairs))
print("Sample of complete pairs:")
print(head(complete_pairs))

# Define all cell types with their proper labels
cell_type_labels <- c(
  # B cells
  "B.cells.naive" = "B cells (naive)",
  "B.cells.memory" = "B cells (memory)",
  "Plasma.cells" = "Plasma cells",
  
  # T cells
  "T.cells.CD8" = "CD8+ T cells",
  "T.cells.CD4.naive" = "CD4+ naive T cells",
  "T.cells.CD4.memory.resting" = "CD4+ memory resting T cells",
  "T.cells.CD4.memory.activated" = "CD4+ memory activated T cells",
  "T.cells.follicular.helper" = "T follicular helper cells",
  "T.cells.regulatory..Tregs." = "Regulatory T cells (Tregs)",
  "T.cells.gamma.delta" = "Gamma T cells",
  
  # NK cells
  "NK.cells.resting" = "NK cells (resting)",
  "NK.cells.activated" = "NK cells (activated)",
  
  # Myeloid cells
  "Monocytes" = "Monocytes",
  "Macrophages.M0" = "Macrophages M0",
  "Macrophages.M1" = "Macrophages M1",
  "Macrophages.M2" = "Macrophages M2",
  "Dendritic.cells.resting" = "Dendritic cells (resting)",
  "Dendritic.cells.activated" = "Dendritic cells (activated)",
  
  # Other immune cells
  "Mast.cells.resting" = "Mast cells (resting)",
  "Mast.cells.activated" = "Mast cells (activated)",
  "Eosinophils" = "Eosinophils",
  "Neutrophils" = "Neutrophils"
)

cell_types <- names(cell_type_labels)

# Define color palette for pre/post treatment
treatment_colors <- c(
  "Pre" = "#2A7D8C",  # Deep Teal
  "Post" = "#B22234"  # Crimson Red
)

# Function to create paired jitter plot
create_paired_jitter_plot <- function(data, cell_type, cell_label = NULL) {
  if (is.null(cell_label)) {
    cell_label <- cell_type
  }
  
  p <- ggplot(data, aes(x = Treatment, y = .data[[cell_type]], group = Patient)) +
    # Add lines connecting paired samples
    geom_line(alpha = 0.5) +
    # Add points
    geom_point(aes(fill = Treatment), 
               size = 4, 
               alpha = 0.7,
               shape = 21,
               stroke = 0.5) +
    # Add styling
    scale_fill_manual(values = treatment_colors) +
    labs(title = paste("Pre vs Post Treatment:", cell_label),
         y = "Cell Proportion",
         x = "Treatment") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}

# Generate paired plots for each cell type
for (cell_type in cell_types) {
  # Create paired plot using complete pairs only
  p <- create_paired_jitter_plot(
    data = complete_pairs,
    cell_type = cell_type,
    cell_label = cell_type_labels[cell_type]
  )
  
  # Save plot
  ggsave(
    filename = file.path(plot_dir, paste0(cell_type, "_paired_jitter.pdf")),
    plot = p,
    width = 6,
    height = 5,
    units = "in",
    dpi = 300
  )
}