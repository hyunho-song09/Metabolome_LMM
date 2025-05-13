
# Load required libraries
library(dplyr)
library(ggplot2)
library(xlsx)
library(mixOmics)
library(vegan)

### 01. Load Input Data ----------------------------------------------------

# Load metabolomics data (Excel)
raw.metabolome <- read.csv("example_data01.csv")
metabolome <- raw.metabolome

# Load microbiome data (CSV)
raw.microbiome <- read.csv("example_data02.csv")
microbiome <- raw.microbiome

### 02. Feature Scaling ----------------------------------------------------

# ----- Metabolome preprocessing -----

# Rename features: met.feature1, met.feature2, ...
n_met_feats <- ncol(metabolome) - 3
met_feat_names <- paste0("met.feature", seq_len(n_met_feats))
colnames(metabolome)[4:ncol(metabolome)] <- met_feat_names

# Scale features (Z-score)
metabolome[, 4:ncol(metabolome)] <- scale(metabolome[, 4:ncol(metabolome)],
                                          center = TRUE, scale = TRUE)

# Optional: mapping table from original names to new feature names
map.metabolome <- data.frame(
  name = colnames(raw.metabolome)[4:ncol(raw.metabolome)],
  index = met_feat_names
)


# ----- Microbiome preprocessing -----

# Rename features: mb.feature1, mb.feature2, ...
n_mb_feats <- ncol(microbiome) - 3
mb_feat_names <- paste0("mb.feature", seq_len(n_mb_feats))
colnames(microbiome)[4:ncol(microbiome)] <- mb_feat_names

# Remove features with all-zero values
feature_cols <- mb_feat_names
non_zero_cols <- feature_cols[colSums(microbiome[, feature_cols]) != 0]
microbiome <- microbiome[, c(colnames(microbiome)[1:3], non_zero_cols)]

# Scale features (Z-score)
microbiome[, 4:ncol(microbiome)] <- scale(microbiome[, 4:ncol(microbiome)],
                                          center = TRUE, scale = TRUE)

# Optional: mapping table from original names to new feature names
map.microbiome <- data.frame(
  name = colnames(raw.microbiome)[4:ncol(raw.microbiome)],
  index = mb_feat_names
)


### 03. sPLS-DA Analysis Pipeline --------------------------------------------

# Define group colors globally
group.colors <- c("Untreated" = "#4CAF50", "Treated" = "#E53935")

# Generalized sPLS-DA function for both week 0 and week 1
run_splsda_pipeline <- function(data,
                                data_type = c("metabolome", "microbiome"),
                                keepX_list,
                                output_prefix,
                                save_fig = FALSE) {
  data_type <- match.arg(data_type)
  
  results <- list()  # To store results for week 0 and 1
  
  for (wk in c(0, 1)) {
    # Filter data by week
    df <- data %>%
      filter(week == wk) %>%
      na.omit() %>%
      mutate(Group = ifelse(treat == "0", "Untreated", "Treated"))
    
    X <- df[, 4:(ncol(df) - 1)]
    Y <- df$Group
    
    # Step 1: keepX tuning
    set.seed(42)
    tune.res <- tune.splsda(X, Y, ncomp = 2,
                            validation = "Mfold", folds = 5, nrepeat = 10,
                            test.keepX = keepX_list, dist = "max.dist",
                            progressBar = TRUE)
    opt.keepX <- tune.res$choice.keepX
    cat(sprintf("[Week %d - %s] Optimal keepX: %s\n", wk, toupper(data_type), opt.keepX))
    
    # Step 2: model fitting
    set.seed(42)
    model <- splsda(X, Y, ncomp = 2, keepX = opt.keepX)
    
    # Step 3: Visualization
    plot_title <- paste("Optimized sPLS-DA (", toupper(data_type), ", Week ", wk, ")", sep = "")
    
    if (save_fig) {
      pdf(paste0(output_prefix, "_Week", wk, "_sPLSDA_", data_type, ".pdf"), width = 6.5, height = 5.5)
      plotIndiv(model, group = Y, ind.names = FALSE, legend = TRUE,
                ellipse = TRUE, col = group.colors,
                title = plot_title)
      dev.off()
    } else {
      # Print directly to RStudio plot viewer
      print(
        plotIndiv(model, group = Y, ind.names = FALSE, legend = TRUE,
                  ellipse = TRUE, col = group.colors,
                  title = plot_title)
      )
    }
    
    # Step 4: Performance Evaluation
    perf.res <- perf(model, validation = "Mfold", folds = 10,
                     nrepeat = 10, progressBar = FALSE)
    err <- perf.res$error.rate
    ber.df <- data.frame(
      Dataset = paste(toupper(data_type), "(optimized, Week ", wk, ")", sep = ""),
      BER.Comp1 = err$BER[1],
      BER.Comp2 = err$BER[2]
    )
    write.csv(ber.df, paste0(output_prefix, "_Week", wk, "_BER_", data_type, ".csv"), row.names = FALSE)
    
    # Step 5: Export selected features
    var1 <- selectVar(model, comp = 1)
    var2 <- selectVar(model, comp = 2)
    df1 <- data.frame(Variable = rownames(var1$value), Loading.Comp1 = var1$value[, 1])
    df2 <- data.frame(Variable = rownames(var2$value), Loading.Comp2 = var2$value[, 1])
    selected.vars <- full_join(df1, df2, by = "Variable")
    write.csv(selected.vars, paste0(output_prefix, "_Week", wk, "_SelectedVars_", data_type, ".csv"), row.names = FALSE)
    
    # Store result
    results[[paste0("week", wk)]] <- list(model = model,
                                          perf = perf.res,
                                          selected = selected.vars,
                                          keepX = opt.keepX)
  }
  
  return(results)
}

# Run for metabolome
splsda.met.res <- run_splsda_pipeline(
  data = metabolome,
  data_type = "metabolome",
  keepX_list = c(10, 15, 20, 30, 40, 60, 100, 150, 200, ncol(metabolome) - 3),
  output_prefix = "250512",
  save_fig = FALSE  # Change to TRUE to export PDF
)

# Run for microbiome
splsda.mb.res <- run_splsda_pipeline(
  data = microbiome,
  data_type = "microbiome",
  keepX_list = c(10, 15, 20, 30, 100, 200, 350, 400, ncol(microbiome) - 3),
  output_prefix = "250512",
  save_fig = FALSE
)


### 04. Volcano Plot  --------------------------------------------

library(tidyverse)
library(ggrepel)

# Generalized volcano plot function
run_volcano_plot <- function(data, mapping_table = NULL, data_type = c("metabolome", "microbiome"),
                             output_prefix = "volcano", save_fig = FALSE) {
  data_type <- match.arg(data_type)
  results.list <- list()
  
  for (w in c(0, 1)) {
    df <- data %>%
      filter(week == w) %>%
      na.omit() %>%
      mutate(Group = ifelse(treat == "0", "Untreated", "Treated"),
             Group = factor(Group, levels = c("Untreated", "Treated")))
    
    df <- cbind(df[, c(1, ncol(df))], df[, 4:(ncol(df) - 1)])
    feature.names <- colnames(df)[3:ncol(df)]
    
    lm_results <- lapply(feature.names, function(var) {
      fit <- lm(as.formula(paste(var, "~ Group")), data = df)
      coef <- summary(fit)$coefficients
      coef["GroupTreated", c("Estimate", "Pr(>|t|)")]
    })
    
    res_lm <- do.call(rbind, lm_results)
    colnames(res_lm) <- c("coef", "pval")
    res_lm <- as.data.frame(res_lm)
    
    # feature 이름 매핑 (metabolome만 해당)
    if (!is.null(mapping_table)) {
      res_lm$feature <- mapping_table$name
    } else {
      res_lm$feature <- rownames(res_lm)
    }
    
    res_lm$week <- w
    results.list[[as.character(w)]] <- res_lm
  }
  
  # Combine results
  res_all <- bind_rows(results.list)
  
  # Volcano grouping
  res_all <- res_all %>%
    mutate(
      direction = case_when(
        coef > 0.2 & pval < 0.05 ~ "UP_strong",
        coef > 0.2 & pval < 0.1 ~ "UP_weak",
        coef < -0.2 & pval < 0.05 ~ "DOWN_strong",
        coef < -0.2 & pval < 0.1 ~ "DOWN_weak",
        TRUE ~ "NO"
      ),
      category = ifelse(direction == "NO", "NO", paste0(week, "W_", direction)),
      delabel = ifelse(direction != "NO", feature, NA)
    )
  
  # Volcano colors
  volcano.colors <- c(
    "0W_UP_strong" = "#1b9e77", "0W_UP_weak" = "#a6d854",
    "0W_DOWN_strong" = "#984ea3", "0W_DOWN_weak" = "#d4b9da",
    "1W_UP_strong" = "#e31a1c", "1W_UP_weak" = "#fdbf6f",
    "1W_DOWN_strong" = "#1f78b4", "1W_DOWN_weak" = "#a6cee3",
    "NO" = "gray70"
  )
  
  # Plotting
  volcano_plot <- ggplot(res_all, aes(x = coef, y = -log10(pval), color = category, label = delabel)) +
    geom_point(size = 2.5, alpha = 0.9) +
    # geom_text_repel(max.overlaps = 10, size = 3) +  # Uncomment to show labels
    scale_color_manual(values = volcano.colors) +
    geom_vline(xintercept = c(-0.2, 0.2), col = "black", linetype = "dotted", size = 0.8) +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dotted", size = 0.8) +
    geom_hline(yintercept = -log10(0.1), col = "black", linetype = "dotted", size = 0.8) +
    xlab("Group (Treated vs Untreated) Coefficient") +
    ylab("-log10(p-value)") +
    ggtitle(paste("Volcano Plot by Week -", toupper(data_type))) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right")
  
  if (save_fig) {
    pdf(paste0(output_prefix, "_", data_type, ".pdf"), width = 8, height = 6.5)
    print(volcano_plot)
    dev.off()
  } else {
    print(volcano_plot)
  }
  
  return(res_all)
}


# Metabolome volcano plot with mapping and figure export
res.met <- run_volcano_plot(
  data = metabolome,
  mapping_table = map.metabolome,
  data_type = "metabolome",
  output_prefix = "volcano_lm_met",
  save_fig = FALSE  # or FALSE to just show in RStudio
)

# Microbiome volcano plot (no mapping, no label)
res.mb <- run_volcano_plot(
  data = microbiome,
  data_type = "microbiome",
  output_prefix = "volcano_lm_mb",
  save_fig = FALSE
)


### 05. LMM Analysis --------------------------------------------

library(lme4)
library(dplyr)

run_lmm_analysis <- function(data, mapping_table = NULL, data_type = c("metabolome", "microbiome")) {
  data_type <- match.arg(data_type)
  n_features <- ncol(data) - 3
  results <- data.frame()
  excluded <- c()
  
  for (i in 1:n_features) {
    feature <- colnames(data)[i + 3]
    subdata <- data[, c("treat", "week", "patient_id", feature)]
    
    if (sum(complete.cases(subdata)) == 0) {
      message("❌ ", feature, ": no non-NA cases. Skipped.")
      excluded <- c(excluded, feature)
      next
    }
    
    tryCatch({
      model <- lmer(as.formula(paste(feature, "~ treat * week + (1 | patient_id)")), data = data)
      coef_table <- summary(model)$coefficients
      for (effect in rownames(coef_table)) {
        results <- rbind(results, data.frame(
          Features = feature,
          Effect = effect,
          Estimate = coef_table[effect, "Estimate"],
          Std.Error = coef_table[effect, "Std. Error"],
          t.value = coef_table[effect, "t value"],
          P.value = 2 * (1 - pnorm(abs(coef_table[effect, "t value"])))
        ))
      }
    }, error = function(e) {
      message("🚨 Error with ", feature, ": ", e$message)
      excluded <- c(excluded, feature)
    })
  }
  
  if (!is.null(mapping_table)) {
    results <- results %>%
      left_join(mapping_table, by = c("Features" = "index")) %>%
      relocate(name, .before = Features)
  } else {
    results$name <- results$Features
  }
  
  # 유의미한 feature 요약
  sig <- results %>%
    filter(Effect != "(Intercept)") %>%
    group_by(name) %>%
    summarise(min_p = min(P.value, na.rm = TRUE)) %>%
    filter(min_p < 0.05)
  
  results_sig <- results %>% filter(name %in% sig$name)
  return(list(all = results, sig = results_sig))
}


### 06. LMM Visualization (Heatmap) --------------------------------------------

library(ComplexHeatmap)
library(reshape)
library(circlize)
library(grid)

plot_lmm_heatmap <- function(data, results_sig, title = "LMM Result",
                             save_fig = FALSE, output_file = "LMM_heatmap.pdf") {
  
  # 1. log-signed p-value
  pval_data <- results_sig %>%
    filter(Effect != "(Intercept)") %>%
    mutate(
      Significance = case_when(
        P.value < 0.01 ~ "**",
        P.value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      log_signed_p = -log10(P.value) * sign(Estimate)
    )
  
  # 2. z-score 평균 요약
  data.z <- data %>% mutate(Group = paste0("treat", treat))
  group_summary <- data.z %>%
    group_by(Group) %>%
    summarise(across(all_of(unique(pval_data$Features)), ~ mean(.x, na.rm = TRUE))) %>%
    as.data.frame()
  
  # 3. 전치 및 이름 매핑
  group_t <- as.data.frame(t(group_summary[, -1]))
  colnames(group_t) <- group_summary$Group
  group_t$Features <- rownames(group_t)
  z_score <- merge(group_t,
                   results_sig[!duplicated(results_sig$name), c("Features", "name")],
                   by = "Features")
  
  # 4. 계수 및 p-value matrix
  es_matrix <- cast(pval_data, name ~ Effect, value = "Estimate")
  pval_matrix <- cast(pval_data, name ~ Effect, value = "P.value")
  
  # 5. 통합
  final <- z_score %>%
    left_join(es_matrix, by = "name") %>%
    left_join(pval_matrix, by = "name", suffix = c("_coef", "_pval"))
  
  final <- final[, c("name", "treat0", "treat1", grep("_coef|_pval", colnames(final), value = TRUE))]
  rownames(final) <- final$name
  final <- final[, -1]
  
  # 6. Heatmap 설정
  col1 <- colorRamp2(c(min(final[,1:2]), 0, max(final[,1:2])), c("#0000CC", "white", "#990000"))
  col2 <- colorRamp2(c(min(final[,3:5])*1.8, 0, max(final[,3:5])*1.8), c("#3399CC", "white", "#FF6666"))
  
  h1 <- Heatmap(data.frame(final[,1:2]), col = col1, cluster_rows = TRUE, show_row_dend = FALSE,
                row_names_side = "left", column_names_side = NULL, border_gp = gpar(col = "black", lwd = 1.2),
                row_km = 2, width = unit(5.5, "cm"), height = unit(15, "cm"))
  
  h2 <- Heatmap(data.frame(final[,3:5]), col = col2, cluster_rows = FALSE, column_split = c(1,2,3),
                width = unit(3, "cm"), height = unit(15, "cm"), row_km = 2,
                cell_fun = function(j, i, x, y, w, h, fill) {
                  grid.rect(x, y, w, h, gp = gpar(col = NA, fill = NA, lwd = 1))
                  pval <- final[,6:8][i, j]
                  if (pval < 0.01) grid.text("**", x, y, gp = gpar(fontsize = 15))
                  else if (pval <= 0.05) grid.text("*", x, y, gp = gpar(fontsize = 15))
                })
  
  if (save_fig) {
    pdf(output_file, width = 8, height = 10)
    print(h1 + h2)
    dev.off()
  } else {
    print(h1 + h2)
  }
}


# Metabolome LMM
lmm.met <- run_lmm_analysis(metabolome, map.metabolome, "metabolome")
plot_lmm_heatmap(metabolome, lmm.met$sig, output_file = "LMM_metabolome.pdf", save_fig = FALSE)

# Microbiome LMM
lmm.mb <- run_lmm_analysis(microbiome, map.microbiome, "microbiome")
plot_lmm_heatmap(microbiome, lmm.mb$sig, output_file = "LMM_microbiome.pdf", save_fig = FALSE)



### 07. Correlation Analysis (Spearman) --------------------------------------------

library(dplyr)
library(ggplot2)
library(xlsx)
library(psych)
library(reshape2)
library(boot)
library(tidyr)
library(stringr)
library(RColorBrewer)

# 1. Load data ---------------------------------------------------------------
raw.df <- read.csv("example_data03.csv")
data <- raw.df
df_list <- split(data, data$Group)

# Extract group-specific subsets
df.week.1 <- rbind(df_list$P_V4, df_list$T_V4)

# 2. Spearman correlation + bootstrapped CI ---------------------------------

# Compute correlation coefficients and p-values
set.seed(42)
w1.results <- corr.test(df.week.1[,5:16], df.week.1[,17:24], method = "spearman")
w1.coef <- melt(w1.results$r)
w1.pval <- melt(w1.results$p)

# Bootstrapping function for confidence interval
boot_spearman <- function(data, indices, var1, var2) {
  d <- data[indices, ]
  cor(d[[var1]], d[[var2]], method = "spearman")
}

# Compute CI for all pairs
ci_list <- list()
x_vars <- colnames(df.week.1)[5:16]
y_vars <- colnames(df.week.1)[17:24]

for (x in x_vars) {
  for (y in y_vars) {
    boot_out <- boot(
      data = df.week.1,
      statistic = function(data, indices) boot_spearman(data, indices, x, y),
      R = 1000
    )
    ci <- boot.ci(boot_out, type = "perc")
    cor_test <- suppressWarnings(cor.test(df.week.1[[x]], df.week.1[[y]], method = "spearman"))
    ci_list[[paste(x, y, sep = "_vs_")]] <- c(
      spearman.rho = cor_test$estimate,
      lower = ci$percent[4],
      upper = ci$percent[5],
      pvalue = cor_test$p.value
    )
  }
}

# Combine results into a dataframe
ci_df <- do.call(rbind, ci_list) %>% as.data.frame()
ci_df$Pair <- rownames(ci_df)
rownames(ci_df) <- NULL
w1.ci_df <- ci_df %>% dplyr::select(Pair, spearman.rho.rho, lower, upper, pvalue)

# 3. Data reshaping for ggplot ----------------------------------------------

ci_df_long.w1 <- w1.ci_df %>%
  separate(Pair, into = c("x_var", "y_var"), sep = "_vs_") %>%
  mutate(
    sig = ifelse(pvalue <= 0.05, "sig", "ns"),
    plot_group = ifelse(sig == "sig", y_var, paste0("ns_", y_var)),
    unified_group = plot_group
  )

# Define custom shapes for key bacteria
shape_list <- c(
  "Lactobacillus_delbrueckii"      = 16,
  "Lachnoclostridium_edouardi"     = 17,
  "Luoshenia_tenuis"               = 15,
  "Faecalibacterium_duncaniae"     = 18,
  "Klebsiella_varicola"            = 3,
  "Anaerovorax_odorimutans"        = 4,
  "Guopingia_tenuis"               = 8,
  "Bacteroides_fragilis"           = 7
)

# Include both significant and non-significant groups
all_groups <- unique(sort(c(ci_df_long.w1$y_var, paste0("ns_", ci_df_long.w1$y_var))))
shape_vals <- setNames(rep(1:25, length.out = length(all_groups)), all_groups)

for (name in names(shape_list)) {
  shape_vals[name] <- shape_list[name]
  shape_vals[paste0("ns_", name)] <- shape_list[name]
}

# Assign colors
set.seed(42)
present_groups <- sort(unique(ci_df_long.w1$unified_group))
color_palette <- c(
  RColorBrewer::brewer.pal(6, "Dark2"),
  rep("lightgray", 9)
)
color_vals <- setNames(color_palette[1:length(present_groups)], present_groups)

# 4. Plotting -----------------------------------------------------------------

ggplot(ci_df_long.w1, aes(y = spearman.rho.rho, x = x_var)) +
  geom_point(
    aes(color = unified_group, shape = unified_group,
        size = ifelse(sig == "sig", "big", "small")),
    position = position_dodge(width = 0.7)
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, color = unified_group),
    width = 1,
    size = 0.8,
    position = position_dodge(width = 0.7)
  ) +
  scale_color_manual(values = color_vals, name = "Bacteria") +
  scale_shape_manual(values = shape_vals, name = "Bacteria") +
  scale_size_manual(values = c("big" = 8, "small" = 2), guide = "none") +
  scale_y_continuous(limits = c(-1, 1), name = "Spearman's rho (± 95% CI)") +
  labs(x = "Metabolites") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )
