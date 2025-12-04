# SoftHybrid

SoftHybrid is an R package for **missing-value imputation in label-free proteomics**, with a particular focus on **single-cell and low-input datasets**.
Instead of treating all missing values as the same, SoftHybrid explicitly models the **continuous spectrum between MAR and MNAR**. It combines:
- **Random Forest (RF)** – to impute *MAR-like* missing values using global structure  
- **Left-shifted minimal probability (minProb)** – to impute *MNAR-like* missing values as low-abundance signals  

These two components are blended using **sigmoid-based weights** that depend on **protein-wise missing rate** and **abundance**, allowing a **smooth transition** between MAR and MNAR regimes rather than a hard threshold.

SoftHybrid was developed and benchmarked on:
- Three-species titration datasets (Dataset HYE) across multiple doses and acquisition modes  
- A single-cell proteomics dataset (Dataset SCP) with known biological structure
- A normoxia vs hypoxia proteomics dataset (Dataset Condition)

## Installation
SoftHybrid is currently available from GitHub.
```r
# install.packages("devtools")
devtools::install_github("YixinShiProteomics/softHybridImpute")
```
Then load the package:
```r
library(softHybridImpute)
```

## Quick Start
### Prepare Inputs
SoftHybrid needs three matrices with identical dimensions and rownames:
- **df_raw:** Raw intensity matrix containing missing values (NA).
- **df_rf:** RF-imputed version of the same matrix (e.g., generated using missForest).
- **df_minProb:** MinProb-imputed version (e.g., using MSnbase or imputeLCMD).

All three must have:
```r
rows = proteins/features
columns = samples
rownames(df_raw) == rownames(df_rf) == rownames(df_minProb)
```
### Run SoftHybrid (single group)
```r
imp <- soft_hybrid_impute(
  df_rf,
  df_minProb,
  df_raw,
  visualize = TRUE,         # saves a TIFF + RDS weight plot
  group_name = "SCP_Group1" # optional label for filenames
)
```
Output:
- A fully imputed matrix (same dimension as input).
- If visualize = TRUE, SoftHybrid also saves:
  
  - GroupName_SoftHybrid_plot.tiff — a scatter plot showing RF weights
  - GroupName_SoftHybrid_plot.rds — ggplot object for editing

## More Information
A detailed description of the SoftHybrid framework, together with full benchmarking analyses, is available in our preprint:

Shi, Y. *et al.* **SoftHybrid: A Hybrid Imputation Algorithm Optimised for Single-Cell Proteomics Data.** bioRxiv, 2025.

## Contact
For questions or suggestions, please open an issue on the GitHub repository:

Issues: https://github.com/YixinShiProteomics/softHybridImpute/issues

Maintainer: Yixin Shi (yixin.shi@wolfson.ox.ac.uk)



