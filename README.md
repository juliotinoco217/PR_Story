# CIBERSORT Immune Cell Analysis by PR Status

This repository contains R code for analyzing and visualizing immune cell composition in breast cancer samples stratified by Progesterone Receptor (PR) status using CIBERSORT results.

## Features

- Processes CIBERSORT deconvolution results for different PR status groups (Low, Medium, High)
- Generates publication-ready violin plots for all immune cell types
- Creates summary statistics for each cell type across PR categories
- Uses a custom color scheme and professional plotting theme

## Cell Types Analyzed

### B cells and Plasma cells
- B cells (naive)
- B cells (memory)
- Plasma cells

### T cells
- CD8+ T cells
- CD4+ naive T cells
- CD4+ memory resting T cells
- CD4+ memory activated T cells
- T follicular helper cells
- Regulatory T cells (Tregs)
- γδ T cells

### NK cells
- NK cells (resting)
- NK cells (activated)

### Myeloid cells
- Monocytes
- Macrophages (M0, M1, M2)
- Dendritic cells (resting and activated)

### Other immune cells
- Mast cells (resting and activated)
- Eosinophils
- Neutrophils

## Output

- PDF violin plots for each cell type in the `plots` directory
- Summary statistics CSV file with mean, standard deviation, and median values

## Dependencies

- ggplot2
- dplyr
- tidyr
