# Probiotics-Induced Metabolome & Microbiome Normalization Study

This repository contains **modular R scripts** and **sample datasets** for analyzing the metabolic and microbial effects of probiotic treatment in patients with liver disease.  
The pipeline focuses on time-dependent and treatment-specific changes using advanced multivariate and longitudinal models.

---
<br>

## Overview

This project aims to investigate:

- Whether **probiotics** normalize the **metabolomic profile** in liver disease patients.
- How the **gut microbiome** changes across **0 and 12 weeks**, depending on **treatment status (Probiotics vs Placebo)**.
- The analysis utilizes **Linear Mixed Models (LMM)** to account for repeated measures and inter-individual variability via `patient_ID` as a random effect.

---
<br>

## How to Run

To execute the full analysis pipeline, simply run:

```r
source("main.R")
```
<br>

## Pipeline Details
<br>

### 1️⃣ Data Loading
📂 **Script:** ```src/01_load_data.R``` <br>
**Description:**

* Loads metabolomics and microbiome datasets (```.csv``` format).

* Stores as ```raw.metabolome```, ```raw.microbiome```.
<br>

### 2️⃣ Feature Preprocessing
📂 **Script:** ```src/02_preprocessing.R``` <br>
**Description:**

* Renames features (```met.feature1```, ```mb.feature1```, etc.).

* Z-score scaling applied to each feature.

* Filters out all-zero features in microbiome data.

* Generates mapping tables for feature tracking.
<br>

### 3️⃣ sPLS-DA (Sparse Partial Least Squares - Discriminant Analysis)
📂 **Script:** ```src/03_splsda.R``` <br>
**Description:**

* Performs classification analysis on both datasets for **week 0** and **week 12**.

* Optimizes ```keepX``` via repeated M-fold cross-validation.

* Plots individual sample distribution by group.

* Outputs feature importance and balanced error rates (BER).
<br>

### 4️⃣ Volcano Plot (Linear Model per Feature)
📂 **Script:** ```src/04_volcano_plot.R``` <br>
**Description:**

* Performs linear regression: ```feature ~ Group``` (Treated vs Untreated).

* Calculates coefficients and p-values.

* Visualizes changes with volcano plots by week.

* Highlights significantly up/downregulated features.
<br>

### 5️⃣ Linear Mixed Model (LMM)
📂 **Script:** ```src/05_lmm_analysis.R``` <br>
**Description:**

* Fits ```feature ~ treat * week + (1 | patient_ID)``` for each feature.

* Extracts effect estimates and p-values.

* Summarizes significant features across all non-intercept effects.

* Uses ```lme4``` package.
<br>

### 6️⃣ Heatmap Visualization of LMM
📂 **Script:** ```src/06_lmm_heatmap.R``` <br>
**Description:**

* Calculates log-signed p-values.

* Aggregates feature-wise z-scores by group.

* Visualizes results using dual heatmaps for z-score and coefficient matrix.

* Highlights statistically significant terms with asterisks.
<br>

### 7️⃣ Spearman Correlation & Confidence Interval
📂 **Script:** ```src/07_correlation_analysis.R``` <br>
**Description:**

* Computes Spearman’s ρ and bootstrapped 95% CI.

* Groups correlation results by significance.

* Visualizes metabolite-microbiome relationships.

* Highlights important microbial features with custom shapes/colors.
<br>

### Example Use Case
This pipeline was applied in a clinical study evaluating **12-week probiotic treatment** in patients with **liver disease**, focusing on:

* **Time × Treatment interaction effects** using LMM

* **Discriminatory features** using sPLS-DA

* **Metabolite-microbiome cross-domain association**
<br>

### 📂 Directory Structure <br>
```
project/
├── main.R
└── src/
    ├── 01_load_data.R
    ├── 02_preprocessing.R
    ├── 03_splsda.R
    ├── 04_volcano_plot.R
    ├── 05_lmm_analysis.R
    ├── 06_lmm_heatmap.R
    └── 07_correlation_analysis.R
```
<br>

### Dependencies
<br>

```
install.packages(c("dplyr", "ggplot2", "boot", "lme4", "ggrepel", 
                   "reshape2", "psych", "tidyr", "stringr", "RColorBrewer"))
if (!requireNamespace("mixOmics")) install.packages("mixOmics")
if (!requireNamespace("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
```
