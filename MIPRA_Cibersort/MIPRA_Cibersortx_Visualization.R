# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Create output directory
plot_dir <- "MIPRA_Cibersortplots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Define plot theme
single_plot_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.background = element_rect(fill = "white"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

# Define cell type labels
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

# Define treatment colors
treatment_colors <- c(
  "Pre" = "#2A7D8C",  # Deep Teal
  "Post" = "#B22234"  # Crimson Red
)

# Read and process CIBERSORT results
cibersort_results <- read.delim("deconvolution/CibersortxResults_MIPRA_.txt", row.names = 1) %>%
  mutate(
    Sample_ID = rownames(.),
    Treatment = factor(ifelse(grepl("B$", Sample_ID), "Pre", "Post"), levels = c("Pre", "Post")),
    Patient = sub("[BS]$", "", Sample_ID)
  ) %>%
  group_by(Patient) %>%
  filter(n() == 2) %>%  # Keep only complete pairs
  ungroup()

# Function to create paired jitter plot
create_paired_jitter_plot <- function(data, cell_type, cell_label = NULL) {
  if (is.null(cell_label)) cell_label <- cell_type
  
  ggplot(data, aes(x = Treatment, y = .data[[cell_type]], group = Patient)) +
    geom_line(alpha = 0.5) +
    geom_point(aes(fill = Treatment), 
               size = 4, 
               alpha = 0.7,
               shape = 21,
               stroke = 0.5) +
    scale_fill_manual(values = treatment_colors) +
    labs(y = paste(cell_label, "(Cell Fraction)"),
         x = NULL) +
    single_plot_theme
}

# Generate plots for each cell type
cell_types <- names(cell_type_labels)
for (cell_type in cell_types) {
  p <- create_paired_jitter_plot(
    data = cibersort_results,
    cell_type = cell_type,
    cell_label = cell_type_labels[cell_type]
  )
  
  ggsave(
    filename = file.path(plot_dir, paste0(cell_type, "_paired_jitter.pdf")),
    plot = p,
    width = 6,
    height = 5,
    units = "in",
    dpi = 300
  )
}

# Response Analysis
###################

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

# Create response mapping dataframe
response_mapping <- data.frame(
  Patient = names(patient_samples),
  Response = sapply(patient_samples, function(x) x$response)
)

# Add response information to the cibersort results
cibersort_with_response <- cibersort_results %>%
  left_join(response_mapping, by = "Patient")

# Function to create paired jitter plot for response groups
create_response_paired_plot <- function(data, cell_type, response_group, cell_label = NULL, show_legend = FALSE) {
  if (is.null(cell_label)) cell_label <- cell_type
  
  # Filter data for specific response group
  plot_data <- data %>%
    filter(Response == response_group)
  
  ggplot(plot_data, aes(x = Treatment, y = .data[[cell_type]], group = Patient)) +
    geom_line(alpha = 0.5, color = "grey50") +
    geom_point(aes(fill = Treatment), 
               size = 4, 
               alpha = 0.7,
               shape = 21,
               stroke = 0.5) +
    scale_fill_manual(values = treatment_colors) +
    labs(y = paste(cell_label, "(Cell Fraction)"),
         x = NULL,
         title = response_group) +
    single_plot_theme +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = if(show_legend) "right" else "none",
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
      # Add grey rectangle around title
      plot.title.background = element_rect(
        fill = "grey90",
        color = "grey50",
        linewidth = 0.5
      )
    )
}

# Generate separate response plots
for (cell_type in cell_types) {
  # Create plots for each response group
  p_nonresponse <- create_response_paired_plot(
    data = cibersort_with_response,
    cell_type = cell_type,
    response_group = "Non-response",
    cell_label = cell_type_labels[cell_type],
    show_legend = FALSE
  )
  
  p_response <- create_response_paired_plot(
    data = cibersort_with_response,
    cell_type = cell_type,
    response_group = "Response",
    cell_label = cell_type_labels[cell_type],
    show_legend = TRUE
  )
  
  # Combine plots side by side (Non-response first, then Response)
  combined_plot <- p_nonresponse + p_response +
    plot_layout(ncol = 2, widths = c(1, 1.2))  # Make right panel slightly wider to accommodate legend
  
  # Save combined plot
  ggsave(
    filename = file.path(plot_dir, paste0(cell_type, "_response_comparison.pdf")),
    plot = combined_plot,
    width = 12,
    height = 5,
    units = "in",
    dpi = 300
  )
}
