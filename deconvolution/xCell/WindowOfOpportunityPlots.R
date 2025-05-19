# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining plots

# Read the xCell results files
mipra_xcell <- read.csv("deconvolution/xCell/MIPRA_xCell_results.csv", row.names = 1)
telepristone_xcell <- read.csv("deconvolution/xCell/Telepristone_GSE126765_xCell_results.csv", row.names = 1)

# Print column names to check available cell types
print("Available cell types in MIPRA dataset:")
print(colnames(mipra_xcell))

# Add a column to identify the dataset source
mipra_xcell$Dataset <- "MIPRA"
telepristone_xcell$Dataset <- "Telepristone"

## Defining Cell Population Lists ##

# CD4 T cells
cd4_populations <- c(
  "CD4..naive.T.cells",
  "CD4..memory.T.cells",
  "CD4..T.cells",
  "CD4..Tcm",
  "CD4..Tem",
  "Th1.cells",
  "Th2.cells",
  "Tregs"
)

# CD8 T cells
cd8_populations <- c(
  "CD8..T.cells",
  "CD8..Tcm",
  "CD8..Tem",
  "CD8..naive.T.cells"
)

# NK cells
nk_populations <- c(
  "NK.cells",
  "NKT"
)

# B cells
b_populations <- c(
  "B.cells",
  "naive.B.cells",
  "pro.B.cells",
  "Memory.B.cells",
  "Class.switched.memory.B.cells",
  "Plasma.cells"
)

# Myeloid cells
myeloid_populations <- c(
  "Monocytes",
  "Macrophages",
  "Macrophages.M1",
  "Macrophages.M2",
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
  "Mast.cells"
)

# Endothelial cells
endothelial_populations <- c(
  "Endothelial.cells",
  "ly.Endothelial.cells",
  "mv.Endothelial.cells"
)

# Epithelial cells
epithelial_populations <- c(
  "Epithelial.cells"
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

## Creating Cell Type-Specific Dataframes ##

# Function to extract cell populations and add cell type label
create_cell_df <- function(data, cell_populations, cell_type_name) {
  df <- data %>%
    select(all_of(c(cell_populations, "Dataset"))) %>%
    pivot_longer(cols = all_of(cell_populations),
                names_to = "Cell_Population",
                values_to = "Score") %>%
    mutate(Cell_Type = cell_type_name)
  return(df)
}

# Create individual dataframes for each cell type
cd4_df <- rbind(
  create_cell_df(mipra_xcell, cd4_populations, "CD4 T cells"),
  create_cell_df(telepristone_xcell, cd4_populations, "CD4 T cells")
)

cd8_df <- rbind(
  create_cell_df(mipra_xcell, cd8_populations, "CD8 T cells"),
  create_cell_df(telepristone_xcell, cd8_populations, "CD8 T cells")
)

nk_df <- rbind(
  create_cell_df(mipra_xcell, nk_populations, "NK cells"),
  create_cell_df(telepristone_xcell, nk_populations, "NK cells")
)

b_df <- rbind(
  create_cell_df(mipra_xcell, b_populations, "B cells"),
  create_cell_df(telepristone_xcell, b_populations, "B cells")
)

myeloid_df <- rbind(
  create_cell_df(mipra_xcell, myeloid_populations, "Myeloid cells"),
  create_cell_df(telepristone_xcell, myeloid_populations, "Myeloid cells")
)

granulocyte_df <- rbind(
  create_cell_df(mipra_xcell, granulocyte_populations, "Granulocytes"),
  create_cell_df(telepristone_xcell, granulocyte_populations, "Granulocytes")
)

endothelial_df <- rbind(
  create_cell_df(mipra_xcell, endothelial_populations, "Endothelial cells"),
  create_cell_df(telepristone_xcell, endothelial_populations, "Endothelial cells")
)

epithelial_df <- rbind(
  create_cell_df(mipra_xcell, epithelial_populations, "Epithelial cells"),
  create_cell_df(telepristone_xcell, epithelial_populations, "Epithelial cells")
)

stromal_df <- rbind(
  create_cell_df(mipra_xcell, stromal_populations, "Stromal cells"),
  create_cell_df(telepristone_xcell, stromal_populations, "Stromal cells")
)

# Create compiled scores dataframe
compiled_df <- rbind(
  create_cell_df(mipra_xcell, compiled_scores, "Compiled Scores"),
  create_cell_df(telepristone_xcell, compiled_scores, "Compiled Scores")
)

# Combine all dataframes into one master dataframe
all_cells_df <- rbind(
  cd4_df,
  cd8_df,
  nk_df,
  b_df,
  myeloid_df,
  granulocyte_df,
  endothelial_df,
  epithelial_df,
  stromal_df,
  compiled_df
)


## Defining Patient Populations by Trial## 

#MIPRA 
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

#Telepristone 
# Create a list of samples with their clinical information
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

#### MIPRA CD8 Spaghetti Plot ### 

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

create_spaghetti_plot <- function(data, x_var, y_var, group_var, 
                                title = NULL, x_lab = NULL, y_lab = NULL,
                                color_palette = NULL) {
  
  # Base plot
  p <- ggplot(data, aes_string(x = x_var, y = y_var, group = group_var)) +
    # Add individual patient lines
    geom_line(color = "gray70", linewidth = 0.5, alpha = 0.6) +
    # Add points for pre and post
    geom_point(aes_string(color = x_var), size = 3) +
    # Facet by cell type
    facet_wrap(~Cell_Type, scales = "free_y", ncol = 2)
  
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

# Define colors for pre and post (keeping Pre first in legend)
mipra_colors <- c(
  "Pre" = "#2A7D8C",   # Deep Teal
  "Post" = "#B22234"   # Crimson Red
)

# Create a data frame for plotting CD8 T cells
mipra_cd8_data <- data.frame()

# Loop through each patient and get their pre/post samples
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
      TimePoint = factor("Pre", levels = c("Pre", "Post")),  # Ensure Pre comes first
      Patient = patient,
      # Clean up the cell type names for display
      Cell_Type = gsub("\\.\\.", "+ ", Cell_Type) %>%
                 gsub("\\.", " ", .)
    )
  
  # Get CD8 data for post samples
  post_data <- mipra_xcell[post_sample, cd8_populations] %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), 
                names_to = "Cell_Type", 
                values_to = "Proportion") %>%
    mutate(
      TimePoint = factor("Post", levels = c("Pre", "Post")),  # Ensure Pre comes first
      Patient = patient,
      # Clean up the cell type names for display
      Cell_Type = gsub("\\.\\.", "+ ", Cell_Type) %>%
                 gsub("\\.", " ", .)
    )
  
  mipra_cd8_data <- rbind(mipra_cd8_data, pre_data, post_data)
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