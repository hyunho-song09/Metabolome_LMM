
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
  
  sig <- results %>%
    filter(Effect != "(Intercept)") %>%
    group_by(name) %>%
    summarise(min_p = min(P.value, na.rm = TRUE)) %>%
    filter(min_p < 0.05)
  
  results_sig <- results %>% filter(name %in% sig$name)
  return(list(all = results, sig = results_sig))
}

