# Drug Sensitivity Explorer: Real GDSC2 Data

## Overview
An interactive Shiny dashboard for exploring drug sensitivity across cancer types using real data from the Genomics of Drug Sensitivity in Cancer (GDSC2) database. The app enables exploration of 286 drugs across 969 cancer cell lines and 31 cancer types, with interactive filtering by drug, cancer type, and targeted pathway.

🚀 **[Open Live App](https://gokul-selvaraj.shinyapps.io/drug-sensitivity-shiny/)**

## Why This Matters
Understanding differential drug sensitivity across cancer types is central to precision oncology and drug development. The GDSC database is one of the most widely used pharmacogenomics resources in cancer research, used to:

- **Identify cancer-type-specific drug sensitivities** — which tumour types respond best to which drugs
- **Discover predictive biomarkers** — correlate genomic features with drug sensitivity
- **Support drug repurposing** — identify new indications for existing compounds
- **Inform clinical trial design** — select patient populations most likely to respond

## Data Source
- **Database:** Genomics of Drug Sensitivity in Cancer (GDSC2)
- **Release:** GDSC2 Release 8.5 (October 2023)
- **Source:** Sanger Institute — https://www.cancerrxgene.org/
- **Scale:** 286 drugs × 969 cancer cell lines × 31 cancer types
- **Measurements:** 195,278 drug sensitivity profiles
- **Metrics:** LN_IC50, AUC, Z-score

## App Features

### Drug Explorer
- Select any of 286 real drugs and visualise sensitivity across all cancer types
- Boxplots ordered by median sensitivity — identify most and least sensitive cancer types
- IC50 distribution histogram with median labelled
- Drug information panel: target, pathway, number of cell lines tested

### Cancer Analysis
- Select any cancer type and identify its most sensitive drugs
- Choose metric: LN_IC50, AUC, or Z-score
- Adjust number of top drugs displayed (5-30)
- Cross-cancer heatmap showing top drugs vs top cancer types

### Pathway Analysis
- Filter by drug target pathway (EGFR signaling, PI3K/mTOR, RTK signaling, etc.)
- Optionally filter by cancer type
- Top 20 most sensitive drugs shown to avoid overcrowding
- Pathway overview showing number of drugs per pathway

### Data Table
- Fully searchable and filterable table of all 195,278 measurements
- Filter by cancer type, drug, and pathway simultaneously
- Download filtered results as CSV

## Technical Stack
- **Framework:** R Shiny + shinydashboard
- **Data:** Real GDSC2 pharmacogenomics data (Sanger Institute)
- **Visualisation:** ggplot2, viridis colour scales
- **Table:** DT (interactive DataTables)
- **Deployment:** shinyapps.io

## How to Run Locally

### Step 1 — Download GDSC2 data
```r
download.file(
  "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.csv",
  "GDSC2_fitted_dose_response_27Oct23.csv"
)
```

### Step 2 — Install R packages
```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr", "tidyr", "DT"))
```

### Step 3 — Run app
```r
setwd("path/to/drug-sensitivity-shiny")
shiny::runApp()
```

Note: The GDSC2 CSV file (~50MB) is excluded from GitHub via .gitignore. Download it separately using the link above.

## Author
**Gokul Selvaraj, PhD**
GitHub: [GokulSelvaraj-Scientist](https://github.com/GokulSelvaraj-Scientist)
Live App: [gokul-selvaraj.shinyapps.io/drug-sensitivity-shiny](https://gokul-selvaraj.shinyapps.io/drug-sensitivity-shiny/)
