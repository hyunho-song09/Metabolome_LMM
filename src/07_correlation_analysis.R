
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