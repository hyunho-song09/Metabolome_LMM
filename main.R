# Source all modular R scripts
source("src/01_load_data.R")
source("src/02_preprocessing.R")
source("src/03_splsda.R")
source("src/04_volcano_plot.R")
source("src/05_lmm_analysis.R")
source("src/06_lmm_heatmap.R")
source("src/07_correlation_analysis.R")

# Run sPLS-DA
splsda.met.res <- run_splsda_pipeline(metabolome, "metabolome", c(10, 15, 20, 30, 40, 60, 100, 150, 200, ncol(metabolome) - 3), "250512", FALSE)
splsda.mb.res <- run_splsda_pipeline(microbiome, "microbiome", c(10, 15, 20, 30, 100, 200, 350, 400, ncol(microbiome) - 3), "250512", FALSE)

# Run Volcano Plot
res.met <- run_volcano_plot(metabolome, map.metabolome, "metabolome", "volcano_lm_met", FALSE)
res.mb <- run_volcano_plot(microbiome, NULL, "microbiome", "volcano_lm_mb", FALSE)

# Run LMM
lmm.met <- run_lmm_analysis(metabolome, map.metabolome, "metabolome")
lmm.mb <- run_lmm_analysis(microbiome, map.microbiome, "microbiome")

# Plot Heatmap
plot_lmm_heatmap(metabolome, lmm.met$sig, "LMM Metabolome", FALSE)
plot_lmm_heatmap(microbiome, lmm.mb$sig, "LMM Microbiome", FALSE)

# Run Correlation Analysis
run_correlation_analysis()
