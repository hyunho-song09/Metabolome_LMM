
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
