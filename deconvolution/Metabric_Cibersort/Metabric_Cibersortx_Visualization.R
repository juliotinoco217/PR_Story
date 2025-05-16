# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining plots

# Set the output directory for plots
plot_dir <- "METABRIC Plots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Read PR stratification data
pr_data <- read.delim("../PR_stratificationIDs.txt")

# Create three dataframes based on PR levels
PR_Low <- pr_data[pr_data$PR_Quartile %in% c("PR_low", "PR_negative"), ]
PR_Medium <- pr_data[pr_data$PR_Quartile == "PR_medium", ]
PR_High <- pr_data[pr_data$PR_Quartile == "PR_high", ]

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

# Set PR_Category as factor with correct order
cibersort_all$PR_Category <- factor(cibersort_all$PR_Category, 
                                   levels = c("Low", "Medium", "High"))

# Create enhanced publication-ready theme
publication_theme <- theme_minimal() + theme(
  # Title styling
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
  # Axis styling
  axis.title.y = element_text(size = 12, face = "bold"),
  axis.title.x = element_blank(),  # Remove x-axis title
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
  legend.position = "none",  # Remove legend since it's redundant
  # Margins and background
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  panel.background = element_rect(fill = "white"),
  strip.background = element_rect(fill = "white"),
  strip.text = element_text(size = 12, face = "bold")
)

# Define color palette
pr_colors <- c(
  "Low" = "#5C6B39",    # Rich Olive Green
  "Medium" = "#2A7D8C", # Deep Teal
  "High" = "#B22234"    # Crimson Red
)

create_violin_plot <- function(data, x_var, y_var, fill_var = NULL, 
                             title = NULL, y_lab = NULL,
                             color_palette = NULL) {
  
  # Base plot
  p <- ggplot(data, aes_string(x = x_var, y = y_var))
  
  # Add fill aesthetic if provided
  if (!is.null(fill_var)) {
    p <- p + aes_string(fill = fill_var)
  }
  
  # Add violin layer
  p <- p + geom_violin(trim = FALSE, 
                      alpha = 0.8,
                      scale = "width", 
                      adjust = 1.2,
                      width = 0.8)
  
  # Add internal boxplot with more visible parameters
  p <- p + geom_boxplot(width = 0.2,
                       alpha = 1,
                       outlier.shape = NA,
                       fill = "white",
                       color = "black",
                       fatten = 1.2)
  
  # Add custom color palette if provided
  if (!is.null(color_palette)) {
    p <- p + scale_fill_manual(values = color_palette)
  }
  
  # Add labels
  p <- p + labs(title = title,
                y = y_lab)
  
  # Y-axis expansion
  p <- p + scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
  
  # Apply publication theme with additional styling
  p <- p + publication_theme + 
    theme(
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3)
    )
  
  return(p)
}

# Generate violin plots for each cell type using the new function
for (cell_type in cell_types) {
  # Create plot title and get proper cell type label
  plot_title <- paste("Distribution of", cell_type_labels[cell_type])
  
  # Create plot using the new function
  p <- create_violin_plot(
    data = cibersort_all,
    x_var = "PR_Category",
    y_var = cell_type,
    fill_var = "PR_Category",
    title = plot_title,
    y_lab = paste(cell_type_labels[cell_type], "Proportion"),
    color_palette = pr_colors
  )
  
  # Save plot
  ggsave(
    filename = file.path(plot_dir, paste0(cell_type, "_violin.pdf")),
    plot = p,
    width = 8,
    height = 7,
    units = "in",
    dpi = 300
  )
}

# Create a summary statistics table
summary_stats <- cibersort_all %>%
  group_by(PR_Category) %>%
  summarise(across(all_of(cell_types), 
                  list(
                    mean = ~mean(., na.rm = TRUE),
                    sd = ~sd(., na.rm = TRUE),
                    median = ~median(., na.rm = TRUE)
                  ))) %>%
  ungroup()

# Save summary statistics
write.csv(summary_stats, file.path(plot_dir, "summary_statistics.csv"), row.names = FALSE)


##### Adaptive Immune Cell Jittered Plot - All subsets in one graph #####

# Step 1: Define adaptive immune cell types and their labels
adaptive_cell_labels <- c(
  # B cells
  "B.cells.naive" = "Naive B cells",
  "B.cells.memory" = "Memory B cells",
  "Plasma.cells" = "Plasma cells",
  
  # T cells
  "T.cells.CD8" = "CD8+ T cells",
  "T.cells.CD4.naive" = "CD4+ naive T cells",
  "T.cells.CD4.memory.resting" = "CD4+ memory\nresting T cells",
  "T.cells.CD4.memory.activated" = "CD4+ memory\nactivated T cells",
  "T.cells.follicular.helper" = "T follicular\nhelper cells",
  "T.cells.regulatory..Tregs." = "Tregs",
  "T.cells.gamma.delta" = "Gamma T cells"
)

adaptive_cell_types <- names(adaptive_cell_labels)

# Step 2: Prepare the data for plotting
adaptive_cell_data <- cibersort_all %>%
  select(PR_Category, all_of(adaptive_cell_types)) %>%
  pivot_longer(
    cols = all_of(adaptive_cell_types),
    names_to = "Cell_Type",
    values_to = "Proportion"
  ) %>%
  mutate(
    Cell_Type = factor(Cell_Type, 
                      levels = adaptive_cell_types,
                      labels = adaptive_cell_labels[adaptive_cell_types])
  )

# Step 3: Create the combined plot
adaptive_cells_plot <- ggplot(adaptive_cell_data, 
                            aes(x = Cell_Type, 
                                y = Proportion, 
                                fill = PR_Category,
                                color = PR_Category)) +
  # Add points with position_jitterdodge
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
             size = 1,
             alpha = 0.5) +
  # Add boxplot layer
  geom_boxplot(position = position_dodge(0.9),
               width = 0.6,
               alpha = 0,  # Transparent fill
               outlier.shape = NA,  # No outlier points (shown in jitter)
               color = "black",
               size = 0.5) +  # Thinner box lines
  # Add styling
  scale_fill_manual(values = pr_colors) +
  scale_color_manual(values = pr_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  labs(title = "Adaptive Immune Cell Distribution by PR Category",
       x = "Cell Type",
       y = "Cell Proportion",
       fill = "PR Category",
       color = "PR Category") +
  single_plot_theme +
  # Additional theme adjustments
  theme(legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# Step 4: Save the plot
ggsave(
  filename = file.path(plot_dir, "adaptive_cells_jitter.pdf"),
  plot = adaptive_cells_plot,
  width = 14,  # Wider to accommodate more cell types
  height = 8,
  units = "in",
  dpi = 300)


##### Innate Immune Cell Jitter Plot - All subsets in one graph #####

# Step 1: Define innate immune cell types and their labels
innate_cell_labels <- c(
  # NK cells
  "NK.cells.resting" = "NK cells\n(resting)",
  "NK.cells.activated" = "NK cells\n(activated)",
  
  # Myeloid cells
  "Monocytes" = "Monocytes",
  "Macrophages.M0" = "M0 Macrophages",
  "Macrophages.M1" = "M1 Macrophages",
  "Macrophages.M2" = "M2 Macrophages",
  "Dendritic.cells.resting" = "Dendritic cells\n(resting)",
  "Dendritic.cells.activated" = "Dendritic cells\n(activated)",
  
  # Other innate cells
  "Mast.cells.resting" = "Mast cells\n(resting)",
  "Mast.cells.activated" = "Mast cells\n(activated)",
  "Eosinophils" = "Eosinophils",
  "Neutrophils" = "Neutrophils"
)

innate_cell_types <- names(innate_cell_labels)

# Step 2: Prepare the data for plotting
innate_cell_data <- cibersort_all %>%
  select(PR_Category, all_of(innate_cell_types)) %>%
  pivot_longer(
    cols = all_of(innate_cell_types),
    names_to = "Cell_Type",
    values_to = "Proportion"
  ) %>%
  mutate(
    Cell_Type = factor(Cell_Type, 
                      levels = innate_cell_types,
                      labels = innate_cell_labels[innate_cell_types])
  )

# Step 3: Create the combined plot
innate_cells_plot <- ggplot(innate_cell_data, 
                           aes(x = Cell_Type, 
                               y = Proportion, 
                               fill = PR_Category,
                               color = PR_Category)) +
  # Add points with position_jitterdodge
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
             size = 1,
             alpha = 0.5) +
  # Add boxplot layer
  geom_boxplot(position = position_dodge(0.9),
               width = 0.6,
               alpha = 0,  # Transparent fill
               outlier.shape = NA,  # No outlier points (shown in jitter)
               color = "black",
               size = 0.5) +  # Thinner box lines
  # Add styling
  scale_fill_manual(values = pr_colors) +
  scale_color_manual(values = pr_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  labs(title = "Innate Immune Cell Distribution by PR Category",
       x = "Cell Type",
       y = "Cell Proportion",
       fill = "PR Category",
       color = "PR Category") +
  single_plot_theme +
  # Additional theme adjustments
  theme(legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# Step 4: Save the plot
ggsave(
  filename = file.path(plot_dir, "innate_cells_jitter.pdf"),
  plot = innate_cells_plot,
  width = 14,  # Slightly wider to accommodate more cell types
  height = 8,
  units = "in",
  dpi = 300)