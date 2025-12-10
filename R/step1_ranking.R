# Step 1: SNP-wise metric calculation and ranking
# ----------------------------------------------
# This function computes several per-SNP metrics using
# univariate logistic regression and returns ranked SNP lists
# for each metric (largest to smallest).

rank_snp_metrics <- function(
    phenotype_train,
    genotype_train,
    phenotype_test,
    genotype_test
) {
  # phenotype_train: numeric vector (0/1), training outcome
  # genotype_train:  matrix/data.frame, rows = individuals, cols = SNPs (training)
  # phenotype_test:  numeric vector (0/1), test/validation outcome
  # genotype_test:   matrix/data.frame, rows = individuals, cols = SNPs (test)
  
  snp_names <- colnames(genotype_train)
  if (is.null(snp_names)) {
    stop("genotype_train must have column names (SNP IDs).")
  }
  
  # Initialize containers
  beta_vec          <- numeric(length(snp_names))
  accuracy_df       <- data.frame(SNP = snp_names, accuracy = 0)
  bal_accuracy_df   <- data.frame(SNP = snp_names, balanced_accuracy = 0)
  f1_df             <- data.frame(SNP = snp_names, f1_score = 0)
  precision_df      <- data.frame(SNP = snp_names, precision = 0)
  pval_df           <- data.frame(SNP = snp_names, minus_log10_p = 0)
  recall_df         <- data.frame(SNP = snp_names, recall = 0)
  roc_auc_df        <- data.frame(SNP = snp_names, roc_auc = 0)
  
  for (i in seq_along(snp_names)) {
    if (i %% 5 == 0) {
      message("Step 1 - processing SNP ", i, " / ", length(snp_names))
    }
    
    df_train <- data.frame(
      y   = phenotype_train,
      snp = genotype_train[, i]
    )
    df_test <- data.frame(
      y   = phenotype_test,
      snp = genotype_test[, i]
    )
    
    # Univariate logistic regression
    mod  <- stats::glm(y ~ snp, data = df_train, family = stats::binomial("logit"))
    summ <- summary(mod)$coefficients
    
    beta_vec[i] <- summ["snp", "Estimate"]
    
    # p-value transformed as -log10(p)
    pval_df$minus_log10_p[i] <- -log10(summ["snp", "Pr(>|z|)"])
    
    # Predicted probabilities on test set for ROC AUC
    prob_test <- stats::predict(mod, newdata = df_test, type = "response")
    pred_obj  <- ROCR::prediction(prob_test, df_test$y)
    auc_obj   <- ROCR::performance(pred_obj, measure = "auc")
    roc_auc_df$roc_auc[i] <- auc_obj@y.values[[1]]
    
    # Directional hard call based on sign of beta
    if (beta_vec[i] >= 0) {
      # 0 = control, 1/2 = case
      y_pred <- ifelse(df_test$snp > 0, 1, 0)
    } else {
      # 0 = case, 1/2 = control
      y_pred <- ifelse(df_test$snp == 0, 1, 0)
    }
    
    # Precision / recall / F1 using ROSE
    measures <- ROSE::accuracy.meas(df_test$y, prob_test, threshold = 0.5)
    precision_df$precision[i] <- measures[[3]]
    recall_df$recall[i]       <- measures[[4]]
    f1_df$f1_score[i]         <- measures[[5]]
    
    # Accuracy and balanced accuracy using yardstick
    df_temp <- data.frame(
      truth     = factor(df_test$y,   levels = c(0, 1)),
      predicted = factor(y_pred,      levels = c(0, 1))
    )
    
    accuracy_df$accuracy[i] <- mean(df_temp$truth == df_temp$predicted)
    
    ba_tbl <- yardstick::bal_accuracy(df_temp, truth = truth, estimate = predicted)
    bal_accuracy_df$balanced_accuracy[i] <- ba_tbl$.estimate[1]
  }
  
  # Rank each metric from largest to smallest
  accuracy_ranked <- accuracy_df[order(-accuracy_df$accuracy), ]
  bal_accuracy_ranked <- bal_accuracy_df[order(-bal_accuracy_df$balanced_accuracy), ]
  f1_ranked <- f1_df[order(-f1_df$f1_score), ]
  precision_ranked <- precision_df[order(-precision_df$precision), ]
  pval_ranked <- pval_df[order(-pval_df$minus_log10_p), ]
  recall_ranked <- recall_df[order(-recall_df$recall), ]
  roc_auc_ranked <- roc_auc_df[order(-roc_auc_df$roc_auc), ]
  
  # Return a named list (names used downstream)
  ranked_list <- list(
    accuracy          = accuracy_ranked,
    balanced_accuracy = bal_accuracy_ranked,
    f1_score          = f1_ranked,
    precision         = precision_ranked,
    pval              = pval_ranked,
    recall            = recall_ranked,
    roc_auc           = roc_auc_ranked
  )
  
  ranked_list
}
