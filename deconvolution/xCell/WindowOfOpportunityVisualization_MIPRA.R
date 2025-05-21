# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining plots
library(readr)  # For read_csv
library(forcats)

# Function to clean column names
clean_colnames <- function(names) {
  # Replace any variation of naive with a standard form
  names <- gsub("na.*ve_B_cells", "naive_B_cells", names, perl = TRUE)
  names <- gsub("naïve_B_cells", "naive_B_cells", names)
  names <- gsub("nave_B_cells", "naive_B_cells", names)
  return(names)
}

# Read and process MIPRA data
mipra_xcell <- read.delim("deconvolution/xCell/MIPRA_xCell_results.txt", 
                         sep = "\t",
                         check.names = FALSE) %>%
  as.data.frame()
rownames(mipra_xcell) <- mipra_xcell$PatientID
mipra_xcell <- mipra_xcell[,-1]  # Remove PatientID column
colnames(mipra_xcell) <- clean_colnames(colnames(mipra_xcell))

# Print column names to check available cell types
print("Available cell types in MIPRA dataset:")
print(colnames(mipra_xcell))

## Define Cell Population Lists ##
cell_populations <- list(
  cd4 = c(
    "CD4_naive_T_cells", "CD4_memory_T_cells", "CD4_T_cells",
    "CD4_Tcm", "CD4_Tem", "Th1_cells", "Th2_cells", "Tregs"
  ),
  cd8 = c(
    "CD8_T_cells", "CD8_Tcm", "CD8_Tem", "CD8_naive_T_cells"
  ),
  nk = c(
    "NK_cells", "NKT"
  ),
  b = c(
    "B_cells", "naive_B_cells", "pro_B_cells", "Memory_B_cells",
    "Class_switched_memory_B_cells", "Plasma_cells"
  ),
  myeloid = c(
    "Monocytes", "Macrophages", "Macrophages_M1", "Macrophages_M2",
    "DC", "cDC", "pDC", "iDC", "aDC"
  ),
  granulocyte = c(
    "Neutrophils", "Basophils", "Eosinophils", "Mast_cells"
  ),
  endothelial = c(
    "Endothelial_cells", "ly_Endothelial_cells", "mv_Endothelial_cells"
  ),
  epithelial = c(
    "Epithelial_cells"
  ),
  stromal = c(
    "Fibroblasts", "MSC", "Pericytes"
  ),
  compiled = c(
    "ImmuneScore", "StromaScore", "MicroenvironmentScore"
  )
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

# Prepare data for plotting
prepare_plot_data <- function(mipra_xcell, patient_samples) {
  tidy_data <- data.frame()
  
  for(patient in names(patient_samples)) {
    pre_sample <- patient_samples[[patient]]$pre
    post_sample <- patient_samples[[patient]]$post
    
    # Get pre and post data
    pre_data <- mipra_xcell[pre_sample, ] %>%
      as.data.frame() %>%
      pivot_longer(cols = everything(), 
                  names_to = "CellType", 
                  values_to = "Score") %>%
      mutate(TimePoint = "Pre", Patient = patient)
    
    post_data <- mipra_xcell[post_sample, ] %>%
      as.data.frame() %>%
      pivot_longer(cols = everything(), 
                  names_to = "CellType", 
                  values_to = "Score") %>%
      mutate(TimePoint = "Post", Patient = patient)
    
    tidy_data <- rbind(tidy_data, pre_data, post_data)
  }
  
  # Calculate delta as Post - Pre for each cell type and patient
  delta <- tidy_data %>%
    pivot_wider(names_from = TimePoint, values_from = Score) %>%
    mutate(Delta = Post - Pre) %>%  # Changed to Post - Pre
    select(Patient, CellType, Delta)
  
  post_data <- tidy_data %>%
    filter(TimePoint == "Post") %>%
    left_join(delta, by = c("Patient", "CellType"))
  
  return(post_data)
}

# Create dot plot
create_category_dotplot <- function(data, category_name, cell_populations) {
  category_data <- data %>%
    filter(CellType %in% cell_populations) %>%
    mutate(
      CellType = gsub("_", " ", CellType),
      Delta_magnitude = abs(Delta)  # Calculate absolute magnitude of change
    )
  
  ggplot(category_data, aes(Patient, fct_reorder(CellType, Delta))) +
    geom_point(aes(size = Delta_magnitude, colour = Delta)) +  # Size now shows magnitude of change
    scale_colour_gradient2(
      high = "#d73027",    # Red for increase (Post > Pre)
      mid = "grey90",      # Grey for no change
      low = "#4575b4",     # Blue for decrease (Post < Pre)
      midpoint = 0
    ) +
    scale_size_continuous(
      range = c(2, 8),
      name = "Magnitude of Change"
    ) +
    labs(
      title = paste(category_name, "Cell Populations"),
      colour = "Change (Post - Pre)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, "lines"),
      legend.position = "right",
      legend.box = "vertical"
    )
}

# Prepare data
plot_data <- prepare_plot_data(mipra_xcell, patient_samples)

# Define plot dimensions
plot_dimensions <- list(
  cd4 = c(10, 6),
  cd8 = c(10, 4),
  nk = c(10, 3),
  b = c(10, 5),
  myeloid = c(10, 7),
  granulocyte = c(10, 4),
  endothelial = c(10, 3),
  epithelial = c(10, 3),
  stromal = c(10, 3),
  compiled = c(10, 3)
)

# Generate and save plots
for(category in names(cell_populations)) {
  plot <- create_category_dotplot(
    plot_data, 
    paste(toupper(substr(category, 1, 1)), substr(category, 2, nchar(category)), sep=""),
    cell_populations[[category]]
  )
  
  ggsave(
    filename = paste0("deconvolution/xCell/plots/MIPRA_", category, "_dotplot.pdf"),
    plot = plot,
    width = plot_dimensions[[category]][1],
    height = plot_dimensions[[category]][2],
    units = "in",
    dpi = 300
  )
}