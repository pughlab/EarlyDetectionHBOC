run_ml_pipeline_repeated <- function(
    data,
    outcome_col,
    features,
    models = select_models,
    folds = 4,
    repeats = 10,
    iterations = 30,
    cohort_name = "cohort",
    outdir_data = ".",
    base_seed = 124
) {
  
  data[[outcome_col]] <- factor(
    data[[outcome_col]],
    levels = c(0, 1),
    labels = c("Negative", "Positive")
  )
  
  all_predictions <- data.frame()
  all_scores      <- data.frame()
  all_importances <- data.frame()
  
  for (iter in 1:iterations) {
    cat("\nIteration:", iter, "\n")
    
    set.seed(base_seed + iter)
    
    train_index <- createDataPartition(
      data[[outcome_col]],
      p = 0.80,
      list = FALSE
    )
    train_data <- data[train_index, ]
    val_data   <- data[-train_index, ]
    
    total_resamples <- folds * repeats
    seeds <- vector("list", total_resamples + 1)
    for (i in seq_len(total_resamples)) {
      seeds[[i]] <- sample.int(1e6, length(models))
    }
    seeds[[total_resamples + 1]] <- sample.int(1e6, 1)
    
    train_control <- trainControl(
      method          = "repeatedcv",
      number          = folds,
      repeats         = repeats,
      classProbs      = TRUE,
      seeds           = seeds,
      summaryFunction = twoClassSummary,
      savePredictions = "final",
      sampling        = "down"
    )
    
    for (model_name in models) {
      cat("\nTraining model:", model_name, "\n")
      formula_str <- paste(outcome_col, "~", paste(features, collapse = " + "))
      
      model <- train(
        as.formula(formula_str),
        data       = train_data,
        method     = model_name,
        preProcess = c("YeoJohnson", "center", "scale"),
        trControl  = train_control,
        metric     = "ROC"
      )
      
      assign(
        paste0("model_", cohort_name, "_", model_name, "_iter_", iter),
        model,
        envir = .GlobalEnv
      )
      
      # Step 1: ROC on training
      train_probs <- predict(model, newdata = train_data, type = "prob")[, "Positive"]
      roc_train <- roc(train_data[[outcome_col]], train_probs)
      
      # Step 2: Get best threshold (Youden's Index)
      best_thresh <- coords(
        roc_train,
        x = "best",
        input = "threshold",
        ret = "threshold",
        best.method = "youden",
        transpose = FALSE
      )[[1]]
      
      # Step 3: Apply to validation set
      val_probs <- predict(model, newdata = val_data, type = "prob")[, "Positive"]
      val_pred  <- factor(
        ifelse(val_probs > best_thresh, "Positive", "Negative"),
        levels = c("Negative", "Positive")
      )
      
      roc_auc <- as.numeric(roc(val_data[[outcome_col]], val_probs)$auc)
      cm      <- confusionMatrix(val_pred, val_data[[outcome_col]])
      
      accuracy    <- cm$overall["Accuracy"]
      precision   <- cm$byClass["Precision"]
      recall      <- cm$byClass["Recall"]
      f1_score    <- 2 * (precision * recall) / (precision + recall)
      kappa       <- cm$overall["Kappa"]
      sensitivity <- cm$byClass["Sensitivity"]
      specificity <- cm$byClass["Specificity"]
      ppv         <- cm$byClass["Pos Pred Value"]
      npv         <- cm$byClass["Neg Pred Value"]
      
      iteration_scores <- data.frame(
        Iteration   = iter,
        Model       = model_name,
        ROC_AUC     = roc_auc,
        Accuracy    = accuracy,
        Precision   = precision,
        Recall      = recall,
        F1_Score    = f1_score,
        Kappa       = kappa,
        Sensitivity = sensitivity,
        Specificity = specificity,
        PPV         = ppv,
        NPV         = npv,
        Threshold   = best_thresh
      )
      
      all_scores      <- bind_rows(all_scores, iteration_scores)
      
      model_predictions <- data.frame(
        model          = model_name,
        iteration      = iter, 
        actual         = val_data[[outcome_col]],
        predicted_prob = val_probs,
        Threshold   = best_thresh
      )
      
      all_predictions <- bind_rows(all_predictions, model_predictions)
      
      if (model_name %in% c("rf", "gbm", "glm")) {
        imp            <- varImp(model)$importance
        imp$Feature    <- rownames(imp)
        imp$Model      <- model_name
        imp$Iteration  <- iter
        all_importances <- bind_rows(all_importances, imp)
      }
    }
  }
  

  # For each model type, find the run whose ROC_AUC is the max

  median_models <- list()
  
  for (model_name in models) {
    df_model_scores <- all_scores[all_scores$Model == model_name, ]
    df_model_scores <- df_model_scores[order(df_model_scores$ROC_AUC), ]
    n_runs          <- nrow(df_model_scores)
    
    median_idx <- which.max(df_model_scores$ROC_AUC)
    median_iter <- df_model_scores$Iteration[median_idx]
    
    var_name    <- paste0(
      "model_",
      cohort_name, "_",
      model_name, "_iter_",
      median_iter
    )
    
    # Also extract the threshold from all_scores
    threshold_for_median <- df_model_scores$Threshold[median_idx]
    
    # Save both model and threshold
    median_models[[model_name]] <- list(
      model     = get(var_name, envir = .GlobalEnv),
      threshold = threshold_for_median
    )
  }
  
  return(list(
    predictions    = all_predictions,
    scores         = all_scores,
    importances    = all_importances,
    median_models  = median_models
  ))
}


convert_results_to_df <- function(scores) {
  scores_summary <- scores %>%
    dplyr::group_by(Model) %>%
    dplyr::summarise(
      N = n(),
      
      Mean_ROC_AUC = mean(ROC_AUC, na.rm = TRUE),
      SD_ROC_AUC = sd(ROC_AUC, na.rm = TRUE),
      Lower_CI_ROC_AUC = Mean_ROC_AUC - 1.96 * SD_ROC_AUC / sqrt(N),
      Upper_CI_ROC_AUC = Mean_ROC_AUC + 1.96 * SD_ROC_AUC / sqrt(N),
      
      Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
      SD_Accuracy = sd(Accuracy, na.rm = TRUE),
      Lower_CI_Accuracy = Mean_Accuracy - 1.96 * SD_Accuracy / sqrt(N),
      Upper_CI_Accuracy = Mean_Accuracy + 1.96 * SD_Accuracy / sqrt(N),
      
      Mean_Precision = mean(Precision, na.rm = TRUE),
      SD_Precision = sd(Precision, na.rm = TRUE),
      Lower_CI_Precision = Mean_Precision - 1.96 * SD_Precision / sqrt(N),
      Upper_CI_Precision = Mean_Precision + 1.96 * SD_Precision / sqrt(N),
      
      Mean_Recall = mean(Recall, na.rm = TRUE),
      SD_Recall = sd(Recall, na.rm = TRUE),
      Lower_CI_Recall = Mean_Recall - 1.96 * SD_Recall / sqrt(N),
      Upper_CI_Recall = Mean_Recall + 1.96 * SD_Recall / sqrt(N),
      
      Mean_F1_Score = mean(F1_Score, na.rm = TRUE),
      SD_F1_Score = sd(F1_Score, na.rm = TRUE),
      Lower_CI_F1_Score = Mean_F1_Score - 1.96 * SD_F1_Score / sqrt(N),
      Upper_CI_F1_Score = Mean_F1_Score + 1.96 * SD_F1_Score / sqrt(N),
      
      Mean_Kappa = mean(Kappa, na.rm = TRUE),
      SD_Kappa = sd(Kappa, na.rm = TRUE),
      Lower_CI_Kappa = Mean_Kappa - 1.96 * SD_Kappa / sqrt(N),
      Upper_CI_Kappa = Mean_Kappa + 1.96 * SD_Kappa / sqrt(N),
      
      Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE),
      SD_Sensitivity = sd(Sensitivity, na.rm = TRUE),
      Lower_CI_Sensitivity = Mean_Sensitivity - 1.96 * SD_Sensitivity / sqrt(N),
      Upper_CI_Sensitivity = Mean_Sensitivity + 1.96 * SD_Sensitivity / sqrt(N),
      
      Mean_Specificity = mean(Specificity, na.rm = TRUE),
      SD_Specificity = sd(Specificity, na.rm = TRUE),
      Lower_CI_Specificity = Mean_Specificity - 1.96 * SD_Specificity / sqrt(N),
      Upper_CI_Specificity = Mean_Specificity + 1.96 * SD_Specificity / sqrt(N),
      
      Mean_FPR = mean(1 - Specificity, na.rm = TRUE),
      SD_FPR = sd(1 - Specificity, na.rm = TRUE),
      Lower_CI_FPR = Mean_FPR - 1.96 * SD_FPR / sqrt(N),
      Upper_CI_FPR = Mean_FPR + 1.96 * SD_FPR / sqrt(N),
      
      Mean_PPV = mean(PPV, na.rm = TRUE),
      SD_PPV = sd(PPV, na.rm = TRUE),
      Lower_CI_PPV = Mean_PPV - 1.96 * SD_PPV / sqrt(N),
      Upper_CI_PPV = Mean_PPV + 1.96 * SD_PPV / sqrt(N),
      
      Mean_NPV = mean(NPV, na.rm = TRUE),
      SD_NPV = sd(NPV, na.rm = TRUE),
      Lower_CI_NPV = Mean_NPV - 1.96 * SD_NPV / sqrt(N),
      Upper_CI_NPV = Mean_NPV + 1.96 * SD_NPV / sqrt(N), 
      
      # === Medians (added below) ===
      Median_PPV = median(PPV, na.rm = TRUE),
      Median_NPV = median(NPV, na.rm = TRUE),
      Median_TPR = median(Sensitivity, na.rm = TRUE),                     # TPR = Sensitivity
      Median_FPR = median(1 - Specificity, na.rm = TRUE)                 # FPR = 1 - Specificity
      
    ) %>%
    dplyr::arrange(desc(Mean_ROC_AUC)) %>%
    dplyr::select(-N)  # Remove if you don’t want sample size in output
  
  return(scores_summary)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plot_roc_curve <- function(results, title, scores_summary_df) {
  results$model <- as.factor(results$model)
  model_levels <- levels(results$model)
  model_colors <- scales::hue_pal()(length(model_levels))
  names(model_colors) <- model_levels
  
  roc_plot <- ggplot()
  
  for (model_name in model_levels) {
    model_data <- results %>% dplyr::filter(model == model_name)
    
    # Compute ROC
    roc_obj <- pROC::roc(
      response = as.numeric(model_data$actual) - 1,
      predictor = model_data$predicted_prob
    )
    
    # Create ROC dataframe
    roc_df <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      model = model_name
    )
    
    # Get AUC from summary dataframe
    auc_summary <- scores_summary_df %>% dplyr::filter(Model == model_name)
    auc_text <- paste0(
      model_name, ": AUC=", round(auc_summary$Mean_ROC_AUC, 3),
      " (95% CI: ", round(auc_summary$Lower_CI_ROC_AUC, 3), "-",
      round(auc_summary$Upper_CI_ROC_AUC, 3), ")"
    )
    
    # Add ROC line
    roc_plot <- roc_plot + 
      geom_line(data = roc_df, aes(x = fpr, y = tpr), color = model_colors[model_name], size = 1.2)
    
    # Add AUC annotation with matching color
    roc_plot <- roc_plot +
      annotate("text",
               x = 0.38,
               y = 0.4 - which(model_levels == model_name) * 0.08,
               label = auc_text,
               size = 4.7,
               hjust = 0,
               color = model_colors[model_name])
  }
  
  # Finalize plot
  roc_plot +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(title = title,
         x = "False Positive Rate",
         y = "True Positive Rate") +
    theme_minimal(base_size= 11) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 14),
      plot.title = element_text(color = "black", face = "bold", size = 13),
      legend.position = "none"
    )
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

evaluate_median_models <- function(median_models, test_data, outcome_col,
                                   outdir_figures = ".", outdir_data = ".",
                                   positive_label = "Positive", negative_label = "Negative") {
  results <- list()
  all_metrics <- data.frame()
  
  for (model_name in names(median_models)) {
    model_obj  <- median_models[[model_name]]$model
    thresh     <- median_models[[model_name]]$threshold
    
    probs <- predict(model_obj, newdata = test_data, type = "prob")[, positive_label]
    true_factor <- factor(
      test_data[[outcome_col]],
      levels = c(0, 1),
      labels = c(negative_label, positive_label)
    )
    
    roc_obj <- pROC::roc(response = true_factor, predictor = probs, levels = c(negative_label, positive_label))
    auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
    
    pred_labels <- factor(
      ifelse(probs > thresh, positive_label, negative_label),
      levels = c(negative_label, positive_label)
    )
    
    cm <- confusionMatrix(pred_labels, true_factor, positive = positive_label)
    
    # Build data.frame with all values inline
    metrics_df <- data.frame(
      Model             = model_name,
      Threshold         = thresh,
      TP                = cm$table["Positive", "Positive"],
      FP                = cm$table["Positive", "Negative"],
      FN                = cm$table["Negative", "Positive"],
      TN                = cm$table["Negative", "Negative"],
      Sensitivity       = round(100 * cm$byClass["Sensitivity"], 2),
      Specificity       = round(100 * cm$byClass["Specificity"], 2),
      PPV               = round(100 * cm$byClass["Pos Pred Value"], 2),
      NPV               = round(100 * cm$byClass["Neg Pred Value"], 2),
      Accuracy          = round(100 * cm$overall["Accuracy"], 2),
      Kappa             = round(cm$overall["Kappa"], 3),
      Balanced_Accuracy = round(100 * cm$byClass["Balanced Accuracy"], 2),
      AUC               = auc_val
    )
    
    all_metrics <- rbind(all_metrics, metrics_df)
    
    # Plot
    plot_cm <- plot_confusion_matrix_theme(
      cm = cm,
      title_label = paste0("", model_name)
    )
    
    # Save plot
    ggsave(file.path(outdir_figures, paste0(model_name, "_confusion_matrix.pdf")),
           plot = plot_cm, width = 6, height = 5)
    
    results[[model_name]] <- list(
      predicted_prob = probs,
      threshold      = thresh,
      confusion      = cm,
      plot           = plot_cm,
      metrics        = metrics_df
    )
  }
  
  # Save one combined CSV
  write.csv(all_metrics,
            file = file.path(outdir_data, "01_classifier_all_models_validation_set_confusion_metrics.csv"),
            row.names = FALSE)
  
  return(results)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot_confusion_matrix_theme <- function(cm, title_label) {
  # Pull out the 2×2 table (rows = Reference/Clinical, cols = Prediction/Plasma)
  tbl <- as.table(cm$table)
  TP <- tbl["Positive","Positive"]
  FP <- tbl["Positive","Negative"]   # Prediction = Positive, Reference = Negative
  FN <- tbl["Negative","Positive"]   # Prediction = Negative, Reference = Positive
  TN <- tbl["Negative","Negative"]
  
  # Build data.frame for ggplot
  plot_df <- data.frame(
    Clinical = factor(
      c("Positive","Negative","Positive","Negative"),
      levels = c("Positive","Negative")
    ),
    Plasma = factor(
      c("Positive","Positive","Negative","Negative"),
      levels = c("Positive","Negative")
    ),
    label = c("TP","FP","FN","TN"),
    Y     = c(TP, FP, FN, TN),
    color = c("good","bad","bad","good")
  )
  
  # Totals by actual (Clinical) class
  row_totals <- c(
    TP + FN,  # total actual Positives
    FP + TN   # total actual Negatives
  )
  
  # Percent of each actual class in each cell
  plot_df$percent <- c(
    round(100 * plot_df$Y[1] / row_totals[1], 2),  # TP/(TP+FN)
    round(100 * plot_df$Y[2] / row_totals[1], 2),  # FN/(TP+FN)
    round(100 * plot_df$Y[3] / row_totals[2], 2),  # FP/(FP+TN)
    round(100 * plot_df$Y[4] / row_totals[2], 2)   # TN/(FP+TN)
  )
  
  plot_df$label_text <- paste0(
    dplyr::recode(plot_df$label,
                  TP = "True positive",
                  FP = "False positive",
                  FN = "False negative",
                  TN = "True negative"),
    "\n(", plot_df$Y, ")\n", plot_df$percent, "%"
  )
  
  # Compute the global metrics
  Sensitivity <- round(100 * TP / (TP + FN), 2)
  Specificity <- round(100 * TN / (TN + FP), 2)
  PPV         <- round(100 * TP / (TP + FP), 2)
  NPV         <- round(100 * TN / (TN + FN), 2)
  
  # Plot
  ggplot(plot_df, aes(x = Clinical, y = Plasma)) +
    geom_tile(
      aes(fill = color),
      color = "black", size = 1.2, alpha = 0.9
    ) +
    geom_text(aes(label = label_text), size = 3, fontface = "bold") +
    scale_fill_manual(values = c("good" = "#b3cde3", "bad" = "#fbb4ae")) +
    xlab(paste0(
      "Clinical Finding\n\n",
      "Sensitivity = ", Sensitivity, "%\n",
      "Specificity = ", Specificity, "%\n",
      "PPV = ", PPV, "% | NPV = ", NPV, "%"
    )) +
    ylab("Plasma Finding") +
    ggtitle(paste0(title_label, " (Samples)")) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(size = 13, hjust = 0.5, face = "bold"),
      axis.title.x    = element_text(face = "bold"),
      axis.title.y    = element_text(face = "bold"),
      axis.text       = element_text(color = "black"),
      legend.position = "none"
    )
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

combine_confusion_plots <- function(results_list,
                                    outdir_figures = ".",
                                    filename = "combined_confusion_matrix.pdf",
                                    ncol = 2,
                                    nrow = 2,
                                    width = 9,
                                    height = 6) {
  
  # Extract plot objects
  plot_list <- lapply(results_list, function(res) res$plot)
  
  # Create grid
  combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol, nrow = nrow)
  
  # Save to file
  output_path <- file.path(outdir_figures, filename)
  ggsave(output_path, plot = combined_plot, width = width, height = height)
  
  return(combined_plot)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot_roc_curves_from_results <- function(results_list,
                                         test_data,
                                         outcome_col,
                                         outdir_figures = ".",
                                         filename_prefix = "model") {
  library(pROC)
  library(ggplot2)
  
  for (model_name in names(results_list)) {
    probs <- results_list[[model_name]]$predicted_prob
    
    true_labels <- factor(test_data[[outcome_col]],
                          levels = c(0, 1),
                          labels = c("Negative", "Positive"))
    
    # Compute ROC and AUC with CI
    roc_obj <- pROC::roc(response = true_labels,
                         predictor = probs,
                         levels = c("Negative", "Positive"))
    
    auc_val <- round(pROC::auc(roc_obj), 3)
    ci_vals <- round(pROC::ci.auc(roc_obj), 3)
    
    auc_text <- paste0("AUC = ", auc_val,
                       " (95% CI: ", ci_vals[1], "-", ci_vals[3], ")")
    
    roc_df <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities
    )
    
    # Plot styled to match your original
    roc_plot <- ggplot(data = roc_df, aes(x = fpr, y = tpr)) +
      geom_line(color = "blue", size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      labs(
        title = paste("ROC Curve -", model_name),
        subtitle = auc_text,
        x = "False Positive Rate",
        y = "True Positive Rate"
      ) +
      theme_minimal() +
      theme(
        panel.background   = element_rect(fill = "white", color = NA),
        plot.background    = element_rect(fill = "white", color = NA),
        text               = element_text(color = "black"),
        axis.text          = element_text(color = "black"),
        axis.title         = element_text(color = "black"),
        plot.title         = element_text(color = "black", size = 10, face = "bold"),
        plot.subtitle      = element_text(color = "black", size = 9),
        legend.background  = element_rect(fill = "white"),
        legend.text        = element_text(color = "black")
      )
    
    # Save
    ggsave(
      filename = file.path(outdir_figures, paste0(filename_prefix, "_", model_name, "_ROC.pdf")),
      plot = roc_plot,
      width = 6,
      height = 5
    )
  }
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Example: print average feature importance for each model
print_feature_importance <- function(importances_df) {
  importances_df %>%
    dplyr::group_by(Model, Feature) %>%
    dplyr::summarise(MeanImportance = mean(Overall, na.rm = TRUE)) %>%
    dplyr::arrange(Model, desc(MeanImportance)) %>%
    print(n = 20)  # change n to see more or fewer rows
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

save_individual_roc_curves <- function(predictions_df, filename_prefix, title, scores_summary_df) {
  predictions_df$model <- as.factor(predictions_df$model)
  
  for (model_name in unique(predictions_df$model)) {
    model_data <- predictions_df %>% dplyr::filter(model == model_name)
    
    # Compute ROC Curve for plotting only
    roc_obj <- pROC::roc(
      response = as.numeric(model_data$actual) - 1,
      predictor = model_data$predicted_prob
    )
    
    roc_df <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities
    )
    
    # Get AUC and 95% CI from summary dataframe
    auc_summary <- scores_summary_df %>% dplyr::filter(Model == model_name)
    auc_text <- paste0("AUC=", round(auc_summary$Mean_ROC_AUC, 3),
                       " (95% CI: ", round(auc_summary$Lower_CI_ROC_AUC, 3), "-",
                       round(auc_summary$Upper_CI_ROC_AUC, 3), ")")
    
    # Generate ROC plot
    roc_plot <- ggplot(data = roc_df, aes(x = fpr, y = tpr)) +
      geom_line(color = "blue", size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      labs(title = paste(title, "-", model_name),
           subtitle = auc_text,
           x = "False Positive Rate",
           y = "True Positive Rate") +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black", size = 10, face = "bold"),
        plot.subtitle = element_text(color = "black", size = 9),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(color = "black")
      )
    
    ggsave(filename = file.path(outdir_figures, paste0(filename_prefix, "_", model_name, "_ROC.pdf")),
           plot = roc_plot, width = 6, height = 5)
  }
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++