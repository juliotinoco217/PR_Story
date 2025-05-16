# Immune Cell Composition Analysis in Breast Cancer by PR Status

## Project Overview
This project investigates the relationship between Progesterone Receptor (PR) status and immune cell composition in breast cancer samples from the METABRIC dataset using CIBERSORTx deconvolution. The analysis examines how immune cell proportions vary across different PR expression levels (Low, Medium, High).

## Repository Structure
```
PR Story/
├── deconvolution/
│   ├── Metabric_Cibersort/
│   │   ├── Metabric_Cibersortx_Visualization.R    # Visualization script
│   │   ├── metabric_cibersortx_statistics.R       # Statistical analysis script
│   │   ├── METABRIC Statistics/                   # Statistical results
│   │   └── METABRIC Plots/                        # Visualization outputs
│   ├── PR_stratificationIDs.txt                   # Sample PR categorization
│   ├── CibersortResults_PR_Low.txt               # CIBERSORT results for PR-low
│   ├── CibersortResults_PRmid.txt                # CIBERSORT results for PR-medium
│   └── CIBERSORTx_Results_PR_high.txt            # CIBERSORT results for PR-high
```

## Methods

### Data Processing
1. **PR Status Stratification**
   - Samples categorized into Low, Medium, and High PR expression groups
   - Stratification based on PR expression quartiles

2. **Immune Cell Deconvolution**
   - Used CIBERSORTx algorithm
   - Analysis of 22 immune cell types:
     * B cells (naive, memory, plasma cells)
     * T cells (CD8+, CD4+ naive/memory, Tfh, Tregs, γδ T cells)
     * NK cells (resting, activated)
     * Myeloid cells (monocytes, M0/M1/M2 macrophages, dendritic cells)
     * Other (mast cells, eosinophils, neutrophils)

### Statistical Analysis
1. **Normality Testing**
   - Shapiro-Wilk test for each cell type and PR category
   - Q-Q plots for visual assessment

2. **Differential Analysis**
   - Kruskal-Wallis test across PR categories
   - Dunn's test for all pairwise comparisons (Low-Medium, Low-High, Medium-High)
   - Benjamini-Hochberg correction for multiple testing
   - Complete analysis for all 22 cell types, regardless of Kruskal-Wallis significance

## Key Findings

### Significant Cell Type Differences (FDR < 0.05)

1. **Most Significant Changes**
   - Macrophages M1 (p_adj = 1.70e-05)
     * Significant decrease in High PR vs both Low and Medium
     * Strongest pairwise difference: Medium-High (statistic = -4.97)
   - Macrophages M2 (p_adj = 4.41e-04)
     * Increased in High PR vs Low (statistic = 4.49)
     * Moderate increase in High vs Medium
   - T cells follicular helper (p_adj = 7.79e-04)
     * Increased in High PR vs Low (statistic = 4.28)
     * Moderate increase in High vs Medium

2. **Other Notable Changes**
   - Plasma cells (p_adj = 2.36e-03)
     * Decreased in High PR vs Low (statistic = -3.93)
   - T cells CD8 (p_adj = 1.11e-02)
     * Decreased in High PR vs Low (statistic = -3.46)
   - T cells CD4 memory activated (p_adj = 2.23e-02)
     * Decreased in High PR vs Low (statistic = -3.08)

### Pattern Analysis
- Most significant differences observed between Low and High PR categories
- Medium PR samples often show intermediate values
- Consistent pattern of altered immune composition with increasing PR expression

## Output Files

### Statistical Results
- `kruskal_wallis_results.csv`: Complete results for all 22 cell types
  * KW statistic, p-value, and adjusted p-value
- `dunns_test_results.csv`: All pairwise comparisons for all cell types
  * Includes all comparisons regardless of significance
  * Statistics, raw p-values, and adjusted p-values
- `significant_differences.csv`: Filtered significant findings (FDR < 0.05)
- `session_info.txt`: R session information for reproducibility

### Visualizations
- Individual boxplots for each cell type showing:
  * Distribution across PR categories
  * Individual data points
  * Kruskal-Wallis test p-value
- Color scheme:
  * Low PR: #5C6B39 (olive)
  * Medium PR: #2A7D8C (teal)
  * High PR: #B22234 (red)

## Dependencies
- R version 4.0 or higher
- Required R packages:
  * ggplot2: Advanced data visualization
  * dplyr: Data manipulation
  * tidyr: Data tidying
  * rstatix: Statistical analysis
  * ggpubr: Publication-ready plots

## Usage

### Running Statistical Analysis
```R
# From the Metabric_Cibersort directory
Rscript metabric_cibersortx_statistics.R
```

### Running Visualizations
```R
# From the Metabric_Cibersort directory
Rscript Metabric_Cibersortx_Visualization.R
```

## Notes
- P-values are adjusted using Benjamini-Hochberg method
- All boxplots include individual data points for transparency
- Non-parametric tests used due to non-normal distribution
- Pairwise comparisons performed for all cell types, providing complete analysis
- Results suggest PR status significantly influences immune composition

## Contact
For questions or issues, please open a GitHub issue or contact the repository maintainer.
