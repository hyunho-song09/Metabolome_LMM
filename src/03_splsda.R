source("src/02_preprocessing.R")
library(mixOmics)

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