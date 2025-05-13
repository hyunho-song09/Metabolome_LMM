# Probiotics-Induced Metabolome & Microbiome Normalization Study

This repository contains **modular R scripts** and **sample datasets** for analyzing the metabolic and microbial effects of probiotic treatment in patients with liver disease.  
The pipeline focuses on time-dependent and treatment-specific changes using advanced multivariate and longitudinal models.

---

## Overview

This project aims to investigate:

- Whether **probiotics** normalize the **metabolomic profile** in liver disease patients.
- How the **gut microbiome** changes across **0 and 12 weeks**, depending on **treatment status (Probiotics vs Placebo)**.
- The analysis utilizes **Linear Mixed Models (LMM)** to account for repeated measures and inter-individual variability via `patient_ID` as a random effect.

---

## How to Run

To execute the full analysis pipeline, simply run:

```r
source("main.R")
