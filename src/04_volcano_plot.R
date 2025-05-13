
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

