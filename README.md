# Drug Sensitivity Prediction: Interactive Shiny Dashboard

## Overview
An interactive R Shiny dashboard for exploring drug sensitivity patterns and predicting treatment response in cancer cell lines. Built on simulated GDSC-style (Genomics of Drug Sensitivity in Cancer) data, this app demonstrates end-to-end analysis from exploratory data visualization through machine learning prediction — all in an interactive interface.

## Background
Predicting which cancer cell lines — and by extension which patients — will respond to a given drug is a central challenge in oncology drug development. The GDSC database is one of the largest public resources for pharmacogenomics data, linking genomic features of cancer cell lines to drug sensitivity measurements (IC50 values). This dashboard simulates that analytical workflow interactively.

## Features

### Overview Tab
- Summary statistics across all 200 cell lines and 8 drugs
- IC50 distribution by cancer type for selected drug
- Drug sensitivity heatmap across all drugs and cancer types
- Mutation frequency profiles by cancer type

### Drug Explorer Tab
- Interactive IC50 distribution with sensitive/resistant classification
- IC50 comparison by mutation status (user-selectable)
- Scatter plots of IC50 vs genomic features (TMB, ploidy, doubling time)
- Summary statistics table by cancer type

### ML Prediction Tab
- Train Random Forest or Logistic Regression models interactively
- Adjustable train/test split
- Live ROC curve with AUC
- Confusion matrix
- Feature importance plot
- Real-time model performance metrics

### Biomarker Analysis Tab
- Statistical testing (Wilcoxon rank-sum) for mutation-drug associations
- Violin plots stratified by cancer type
- Biomarker effect heatmap across all drug-mutation combinations

## Dataset
- **Simulated** GDSC-style data based on published distributions
- 200 cancer cell lines across 8 cancer types
- 8 drugs with realistic genomic associations (e.g., EGFR mutation → Erlotinib sensitivity)
- 12 genomic features: TP53, KRAS, BRAF, EGFR, PIK3CA mutations, MYC amplification, PTEN loss, MSI status, TMB, ploidy, doubling time

## How to Run

### Locally
1. Install R and RStudio
2. Install required packages:
```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr", "tidyr",
                   "randomForest", "caret", "pROC", "RColorBrewer", "DT"))
```
3. Run the app:
```r
shiny::runApp("app.R")
```

### Online
Deploy to shinyapps.io (free tier available):
```r
install.packages("rsconnect")
rsconnect::deployApp(".")
```

## Technical Stack
- **Frontend:** R Shiny + shinydashboard
- **ML:** caret, randomForest, pROC
- **Visualization:** ggplot2, DT
- **Statistics:** Wilcoxon rank-sum test

## Author
**Gokul Selvaraj, PhD**
GitHub: [GokulSelvaraj-Scientist](https://github.com/GokulSelvaraj-Scientist)
