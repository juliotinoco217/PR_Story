# Immune Cell Composition Analysis in Breast Cancer by PR Status

## Project Overview
This project investigates the relationship between Progesterone Receptor (PR) status and immune cell composition in breast cancer samples from the METABRIC dataset using CIBERSORTx deconvolution.

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
   - 22 immune cell types analyzed
   - Results separated by PR status groups

### Statistical Analysis
1. **Normality Testing**
   - Shapiro-Wilk test for each cell type and PR category
   - Q-Q plots for visual assessment

2. **Differential Analysis**
   - Kruskal-Wallis test across PR categories
   - Dunn's test for post-hoc pairwise comparisons
   - Benjamini-Hochberg correction for multiple testing

## Key Findings

### Significant Cell Type Differences (FDR < 0.05)

1. **Most Significant Changes**
   - Macrophages M1 (p_adj = 1.70e-05)
     * Decreased in High PR vs both Low and Medium
   - Macrophages M2 (p_adj = 4.41e-04)
     * Increased in High PR vs Low
   - T cells follicular helper (p_adj = 7.79e-04)
     * Increased in High PR vs Low

2. **Other Significant Changes**
   - Plasma cells (p_adj = 2.36e-03)
     * Decreased in High PR vs Low
   - T cells CD8 (p_adj = 1.11e-02)
     * Decreased in High PR vs Low
   - T cells CD4 memory activated (p_adj = 2.23e-02)
     * Decreased in High PR vs Low
   - Dendritic cells resting (p_adj = 2.23e-02)
     * Different between Low-Medium

### Pairwise Comparisons
Most significant differences were observed between Low and High PR categories, suggesting a gradient effect of PR expression on immune composition.

## Output Files

### Statistical Results
- `kruskal_wallis_results.csv`: Overall differences between groups
- `dunns_test_results.csv`: Pairwise comparison results
- `significant_differences.csv`: Significant findings (FDR < 0.05)
- `session_info.txt`: R session information for reproducibility

### Visualizations
- Individual boxplots for each cell type
- Q-Q plots for normality assessment
- Combined visualization plots

## Dependencies
- R version 4.0 or higher
- Required R packages:
  * ggplot2
  * dplyr
  * tidyr
  * rstatix
  * ggpubr

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
- Non-parametric tests were used due to non-normal distribution of data

## Contact
For questions or issues, please open a GitHub issue or contact the repository maintainer.
