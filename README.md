# Visualizations
R scripts for data visualization using heatmaps, violin plots, and barplots

# Clinical Data Visualization in R

This repository contains scripts to generate:
- Complex heatmaps with layered annotations
- Grouped barplots with statistical tests
- Violin plots (to be added)

## 📦 Requirements

Install the following R packages:

```R
install.packages(c("ggplot2", "ggpubr", "dplyr", "openxlsx"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "circlize"))

