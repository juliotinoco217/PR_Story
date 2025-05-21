# MIPRA CIBERSORT Analysis Pipeline

## Overview
This repository contains an R-based analysis pipeline for processing and analyzing CIBERSORT immune cell deconvolution data from the MIPRA study. The pipeline performs statistical analyses to identify differences in immune cell populations between pre- and post-treatment samples, as well as between response groups.

## Data Structure

### Input Data
- **Location**: `deconvolution/CibersortxResults_MIPRA_.txt`
- **Format**: Tab-delimited text file
- **Content**: CIBERSORT results containing immune cell proportions
- **Required Columns**:
  - Row names: Sample IDs (format: `MXXXB` for pre-treatment, `MXXXS` for post-treatment)
  - Columns: 22 immune cell type proportions
  - Additional columns: P-value, Correlation, RMSE (CIBERSORT metrics)

### Sample Information
The analysis includes 8 paired samples (pre/post treatment):
- **Responders** (n=4):
  - M062 (B/S)
  - M070 (B/S)
  - M077 (B/S)
  - M140 (B/S)
- **Non-responders** (n=4):
  - M073 (B/S)
  - M090 (B/S)
  - M094 (B/S)
  - M105 (B/S)

### Immune Cell Types Analyzed
1. B cell populations:
   - B cells naive
   - B cells memory
   - Plasma cells

2. T cell populations:
   - T cells CD8
   - T cells CD4 naive
   - T cells CD4 memory resting
   - T cells CD4 memory activated
   - T cells follicular helper
   - T cells regulatory (Tregs)
   - T cells gamma delta

3. NK cell populations:
   - NK cells resting
   - NK cells activated

4. Myeloid cell populations:
   - Monocytes
   - Macrophages M0
   - Macrophages M1
   - Macrophages M2
   - Dendritic cells resting
   - Dendritic cells activated

5. Other immune cells:
   - Mast cells resting
   - Mast cells activated
   - Eosinophils
   - Neutrophils

## Analysis Pipeline

### Statistical Methods

1. **Pre vs Post Treatment Analysis**
   - Method: Paired t-test
   - Comparison: Post-treatment vs Pre-treatment for each cell type
   - Metrics reported:
     - t-statistic
     - P-value
     - Mean difference
     - 95% Confidence Interval
     - FDR-adjusted P-value

2. **Response Group Analysis**
   - Method: Independent t-test
   - Comparisons: 
     - Responders vs Non-responders at pre-treatment
     - Responders vs Non-responders at post-treatment
   - Metrics reported:
     - t-statistic
     - P-value
     - Mean difference
     - 95% Confidence Interval
     - FDR-adjusted P-value

### Multiple Testing Correction
- Method: Benjamini-Hochberg False Discovery Rate (FDR)
- Applied separately for:
  - Pre vs Post comparisons
  - Response group comparisons

## Output Files

### Location
All output files are saved in the `MIPRA_Statistics` directory with timestamps.

### File Types
1. **Paired t-test Results** (`paired_ttest_results_[timestamp].csv`):
   - Cell_Type: Name of the immune cell population
   - Statistic: t-statistic
   - P_Value: Unadjusted p-value
   - Mean_Diff: Mean difference (Post - Pre)
   - CI_Lower: Lower bound of 95% CI
   - CI_Upper: Upper bound of 95% CI
   - P_Adjusted: FDR-adjusted p-value

2. **Response Group Results** (`response_ttest_results_[timestamp].csv`):
   - Cell_Type: Name of the immune cell population
   - Timepoint: Pre or Post treatment
   - Statistic: t-statistic
   - P_Value: Unadjusted p-value
   - Mean_Diff: Mean difference (Response - Non-response)
   - CI_Lower: Lower bound of 95% CI
   - CI_Upper: Upper bound of 95% CI
   - P_Adjusted: FDR-adjusted p-value

## Dependencies
Required R packages:
- `ggplot2`: For potential visualization extensions
- `dplyr`: For data manipulation
- `tidyr`: For data reshaping

## Usage
1. Ensure input file is present in `deconvolution/CibersortxResults_MIPRA_.txt`
2. Run the analysis script:
   ```R
   Rscript MIPRA_Cibersort/MIPRA_Cibersort_Stats.R
   ```
3. Results will be saved in `MIPRA_Statistics` directory
4. Significant results (P < 0.05) will be printed to console

## Interpretation Guidelines
1. **P-values**:
   - Nominal significance: P < 0.05
   - Multiple testing corrected significance: FDR < 0.05

2. **Effect Sizes**:
   - Mean differences represent absolute changes in cell proportions
   - Positive values indicate:
     - For paired tests: Higher in post-treatment
     - For response tests: Higher in responders

3. **Confidence Intervals**:
   - Non-overlapping zero indicates significant difference
   - Width indicates precision of the estimate

## Notes and Limitations
1. Small sample size (n=8 pairs) limits statistical power
2. Multiple testing correction may be conservative
3. Analysis assumes:
   - Normal distribution of differences (paired t-test)
   - Equal variances between groups (response analysis)
   - Independence of samples

## Contact
For questions or issues, please open a GitHub issue or contact the repository maintainer.
