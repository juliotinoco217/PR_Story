# PR Status and Immune Cell Composition Analysis

This repository contains code and analysis for investigating the relationship between Progesterone Receptor (PR) status and immune cell composition in breast cancer samples from the METABRIC dataset.

## Project Structure

```
PR Story/
├── deconvolution/
│   └── Metabric_Cibersort/
│       ├── Metabric_Cibersortx_Visualization.R
│       └── METABRIC Plots/
├── PR_stratificationIDs.txt
├── CibersortResults_PR_Low.txt
├── CibersortResults_PRmid.txt
└── CIBERSORTx_Results_PR_high.txt
```

## Analysis Components

### Immune Cell Deconvolution
- Uses CIBERSORTx results to analyze immune cell composition
- Stratifies samples based on PR status (Low/Negative, Medium, High)
- Includes both adaptive and innate immune cell populations

### Visualizations
The R script (`Metabric_Cibersortx_Visualization.R`) generates several plots:

1. **Adaptive Immune Cell Distribution**
   - B cells (naive, memory, and plasma cells)
   - T cells (CD8+, CD4+ naive, memory resting/activated, T follicular helper, Tregs, and Gamma T cells)
   - Shows distribution across PR categories with boxplots and individual data points

2. **Innate Immune Cell Distribution**
   - NK cells (resting and activated)
   - Myeloid cells (Monocytes, M0/M1/M2 Macrophages, Dendritic cells)
   - Other innate cells (Mast cells, Eosinophils, Neutrophils)
   - Displays distribution patterns across PR categories

### Output Files
- Individual cell type plots: `*_violin.pdf`
- Combined adaptive immune cell plot: `adaptive_cells_jitter.pdf`
- Combined innate immune cell plot: `innate_cells_jitter.pdf`
- Summary statistics: `summary_statistics.csv`

## Dependencies
- R packages:
  - ggplot2
  - dplyr
  - tidyr
  - patchwork

## Usage
1. Ensure all required R packages are installed
2. Place CIBERSORTx result files in the appropriate directory
3. Run the R script:
```R
Rscript deconvolution/Metabric_Cibersort/Metabric_Cibersortx_Visualization.R
```
