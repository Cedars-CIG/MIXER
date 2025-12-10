# Step 2: Ridge regression, metric-level weights, SNP-level adaptive weights
# -------------------------------------------------------------------------


# 2a. Ridge regression helper
Ridge_func <- function(y, geno, family = "binomial", nfolds = 20) {
  # y: response vector (e.g., 0/1 for binomial)
  # geno: genotype matrix or data.frame (rows: samples, cols: SNPs)
  # family: glmnet family, default "binomial"
  # nfolds: number of CV folds for cv.glmnet
  
  x <- as.matrix(geno)
  
  message("Step 2a - Ridge regression start...")
  
  type_measure <- if (family == "binomial") "deviance" else "mse"
  
  cvfit <- glmnet::cv.glmnet(
    x      = x,
    y      = y,
    family = family,
    alpha  = 0,              # ridge
    type.measure = type_measure,
    intercept    = TRUE,
    nfolds       = nfolds
  )
  
  fit <- glmnet::glmnet(
    x      = x,
    y      = y,
    family = family,
    lambda = cvfit$lambda.min,
    alpha  = 0,
    intercept = TRUE
  )
  
  message("Step 2a - Ridge regression done")
  
  coefs_mat <- stats::coef(fit, s = cvfit$lambda.min)
  coefs_mat <- as.matrix(coefs_mat)
  
  df_coef <- data.frame(
    SNP  = rownames(coefs_mat),
    beta = as.numeric(coefs_mat),
    row.names = NULL
  )
  
  # remove intercept row
  df_coef <- df_coef[df_coef$SNP != "(Intercept)", , drop = FALSE]
  
  df_coef
}


# 2b. Threshold-moving helper: performance vs threshold on validation set
threshold_moving_func <- function(mod, df, threshold_list = NULL, n_grid = 1000) {
  # mod: fitted glm model, usually y ~ prs, family = binomial
  # df:  data.frame with columns y (0/1) and prs (numeric PRS)
  # threshold_list: optional numeric vector of thresholds;
  #                 if NULL, a grid is generated from prob range.
  # n_grid: number of points if threshold_list is not supplied.
  
  prob <- stats::predict(mod, newdata = df, type = "response")
  
  # AUC on the provided df (validation)
  pred_obj <- ROCR::prediction(prob, df$y)
  auc_val  <- ROCR::performance(pred_obj, measure = "auc")@y.values[[1]]
  
  # p-value for prs (from model)
  summ <- summary(mod)$coefficients
  if (!"prs" %in% rownames(summ)) {
    stop("Model in threshold_moving_func() must include 'prs' as a coefficient.")
  }
  pval <- summ["prs", "Pr(>|z|)"]
  pval <- -log10(pval)
  
  if (is.null(threshold_list)) {
    prob_min <- min(prob, na.rm = TRUE)
    prob_max <- max(prob, na.rm = TRUE)
    threshold_list <- seq(prob_min, prob_max, length.out = n_grid)
  }
  
  message("threshold moving start...")
  
  rslt <- data.frame(
    threshold       = threshold_list,
    above_threshold = NA_real_,
    accuracy        = NA_real_,
    ba              = NA_real_,
    precision       = NA_real_,
    recall          = NA_real_,
    F1              = NA_real_,
    auc             = rep(auc_val, length(threshold_list)),
    pval            = rep(pval,    length(threshold_list))
  )
  
  for (i in seq_along(threshold_list)) {
    if (i %% 500 == 0) {
      message("  threshold ", i, " / ", length(threshold_list), " done.")
    }
    
    thr <- threshold_list[i]
    rslt$above_threshold[i] <- mean(prob >= thr)
    
    # Precision / recall / F1 using ROSE
    measures <- ROSE::accuracy.meas(df$y, prob, threshold = thr)
    rslt$precision[i] <- measures[[3]]
    rslt$recall[i]    <- measures[[4]]
    rslt$F1[i]        <- measures[[5]]
    
    # Hard classification for accuracy and balanced accuracy
    y_pred <- as.integer(prob >= thr)
    if (length(unique(y_pred)) == 2) {
      df_temp <- data.frame(
        truth     = factor(df$y,    levels = c(0, 1)),
        predicted = factor(y_pred,  levels = c(0, 1))
      )
      
      rslt$accuracy[i] <- mean(df_temp$truth == df_temp$predicted)
      ba_tbl           <- yardstick::bal_accuracy(df_temp, truth = truth, estimate = predicted)
      rslt$ba[i]       <- ba_tbl$.estimate[1]
    }
  }
  
  rslt
}


# 2c. Compute metric-level weights from PRS performance (train/validation)
compute_metric_weights <- function(
    y_train,
    y_val,
    geno_train_list,
    geno_val_list,
    ridge_coef_list,
    n_threshold = 2000
) {
  # y_train: outcome vector for training set (0/1)
  # y_val:   outcome vector for validation set (0/1)
  # geno_train_list: *named* list of training genotype matrices (one per metric)
  # geno_val_list:   *named* list of validation genotype matrices (one per metric)
  # ridge_coef_list: *named* list of ridge coef data.frames with columns SNP, beta
  # n_threshold: number of thresholds for scanning in threshold_moving_func
  
  if (is.null(names(geno_train_list)) ||
      is.null(names(geno_val_list))   ||
      is.null(names(ridge_coef_list))) {
    stop("geno_train_list, geno_val_list, and ridge_coef_list must be *named* lists.")
  }
  
  metrics <- names(geno_train_list)
  
  if (!identical(sort(metrics), sort(names(geno_val_list))) ||
      !identical(sort(metrics), sort(names(ridge_coef_list)))) {
    stop("Names of geno_train_list, geno_val_list, and ridge_coef_list must match.")
  }
  
  # Reorder lists to consistent order
  geno_train_list <- geno_train_list[metrics]
  geno_val_list   <- geno_val_list[metrics]
  ridge_coef_list <- ridge_coef_list[metrics]
  
  # 1) For each metric: build PRS, fit glm, run threshold-moving
  rslt_list <- vector("list", length(metrics))
  names(rslt_list) <- metrics
  
  for (i in seq_along(metrics)) {
    crit <- metrics[i]
    message("Step 2b - Processing metric / SNP set: ", crit)
    
    geno_train <- geno_train_list[[i]]
    geno_val   <- geno_val_list[[i]]
    coef_df    <- ridge_coef_list[[i]]  # data.frame(SNP, beta)
    
    # Align SNPs between geno and coef
    snp_common <- intersect(colnames(geno_train), coef_df$SNP)
    
    if (length(snp_common) == 0) {
      warning("No overlapping SNPs for metric ", crit, ". Skipping; result will be NULL.")
      rslt_list[[i]] <- NULL
      next
    }
    
    coef_vec <- coef_df$beta[match(snp_common, coef_df$SNP)]
    
    x_train <- as.matrix(geno_train[, snp_common, drop = FALSE])
    x_val   <- as.matrix(geno_val[,   snp_common, drop = FALSE])
    
    # PRS
    prs_train <- as.numeric(x_train %*% coef_vec)
    prs_val   <- as.numeric(x_val   %*% coef_vec)
    
    df_train <- data.frame(
      y   = y_train,
      prs = prs_train
    )
    df_val <- data.frame(
      y   = y_val,
      prs = prs_val
    )
    
    # Fit glm on training PRS
    mod <- stats::glm(y ~ prs, data = df_train, family = stats::binomial("logit"))
    
    # Threshold grid based on training probabilities
    prob_train     <- stats::predict(mod, newdata = df_train, type = "response")
    prob_min       <- min(prob_train, na.rm = TRUE)
    prob_max       <- max(prob_train, na.rm = TRUE)
    threshold_list <- seq(prob_min, prob_max, length.out = n_threshold)
    
    # Threshold moving on validation
    rslt_list[[i]] <- threshold_moving_func(
      mod,
      df_val,
      threshold_list = threshold_list
    )
  }
  
  # Drop NULL entries
  valid_idx  <- !vapply(rslt_list, is.null, logical(1))
  rslt_list  <- rslt_list[valid_idx]
  metrics    <- names(rslt_list)
  
  if (length(rslt_list) == 0) {
    stop("No valid threshold-moving results; check inputs to compute_metric_weights().")
  }
  
  message("Step 2b - summary start...")
  
  # 2) Summarize rslt_list into metric_weights
  prevalence <- mean(y_val)
  
  best_metric <- data.frame(
    metrics         = metrics,
    threshold       = numeric(length(metrics)),
    above_threshold = numeric(length(metrics)),
    accuracy        = numeric(length(metrics)),
    ba              = numeric(length(metrics)),
    F1              = numeric(length(metrics)),
    precision       = numeric(length(metrics)),
    pval            = numeric(length(metrics)),
    recall          = numeric(length(metrics)),
    auc             = numeric(length(metrics)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(metrics)) {
    crit <- metrics[i]
    df   <- rslt_list[[i]]
    df   <- stats::na.omit(df)
    
    # Keep rows with above_threshold >= prevalence
    df <- df[df$above_threshold >= prevalence, , drop = FALSE]
    if (nrow(df) == 0) {
      warning("No thresholds with above_threshold >= prevalence for ", crit,
              ". Metrics will remain 0 for this set.")
      next
    }
    
    df$above_per_distance <- abs(df$above_threshold - prevalence)
    threshold_best <- df$threshold[
      df$above_per_distance == min(df$above_per_distance, na.rm = TRUE)
    ]
    
    if (length(threshold_best) > 1) {
      threshold_best <- max(threshold_best, na.rm = TRUE)
    }
    
    message(crit)
    message("  threshold_best: ", threshold_best)
    
    best_row <- df[df$threshold == threshold_best, , drop = FALSE]
    
    best_metric$threshold[i]       <- threshold_best
    best_metric$above_threshold[i] <- best_row$above_threshold
    best_metric$accuracy[i]        <- best_row$accuracy
    best_metric$ba[i]              <- best_row$ba
    best_metric$precision[i]       <- best_row$precision
    best_metric$pval[i]            <- best_row$pval
    best_metric$recall[i]          <- best_row$recall
    best_metric$F1[i]              <- best_row$F1
    best_metric$auc[i]             <- best_row$auc
  }
  
  final <- best_metric[, -c(2, 3), drop = FALSE]  # drop threshold + above_threshold
  colnames(final) <- c("SNP_set", "accuracy", "balanced_accuracy",
                       "f1_score", "precision", "pval", "recall", "roc_auc")
  
  # Min–max normalization per metric column
  normalized_metric <- final
  for (j in 2:ncol(normalized_metric)) {
    col     <- normalized_metric[[j]]
    col_min <- min(col, na.rm = TRUE)
    col_max <- max(col, na.rm = TRUE)
    
    if (isTRUE(all.equal(col_min, col_max))) {
      normalized_metric[[j]] <- 0
    } else {
      normalized_metric[[j]] <- (col - col_min) / (col_max - col_min)
    }
  }
  
  # Aggregate across metrics
  M_minmax <- rowSums(normalized_metric[, -1, drop = FALSE], na.rm = TRUE)
  
  metric_weights <- data.frame(
    metric   = normalized_metric$SNP_set,
    M_minmax = M_minmax,
    stringsAsFactors = FALSE
  )
  
  metric_weights
}


# 2d. Combine ridge coefs + metric weights into SNP-level adaptive weights
compute_snp_weights <- function(ridge_coef_list, metric_weights, epsilon = 1e-4) {
  # ridge_coef_list: *named* list of data.frames, each with columns SNP, beta
  # metric_weights:  data.frame with columns:
  #                  - metric   (matching names(ridge_coef_list))
  #                  - M_minmax (combined performance score)
  # epsilon: small constant to avoid zero weights
  
  if (is.null(names(ridge_coef_list))) {
    stop("ridge_coef_list must be a *named* list.")
  }
  
  metrics <- names(ridge_coef_list)
  
  if (!all(metrics %in% metric_weights$metric)) {
    stop("All names(ridge_coef_list) must appear in metric_weights$metric.")
  }
  
  # Union of all SNPs across metrics
  all_snps <- sort(unique(unlist(lapply(ridge_coef_list, function(df) df$SNP))))
  df_coef  <- data.frame(SNP = all_snps, stringsAsFactors = FALSE)
  
  # Fill one column per metric with |beta|
  for (m in metrics) {
    df_m <- ridge_coef_list[[m]]
    
    col_vals <- numeric(length(all_snps))
    idx      <- match(df_m$SNP, all_snps)
    col_vals[idx] <- df_m$beta
    
    df_coef[[m]] <- abs(col_vals)
  }
  
  # Metric weights in same order as metrics
  w <- metric_weights$M_minmax[match(metrics, metric_weights$metric)]
  
  beta_mat     <- as.matrix(df_coef[, metrics, drop = FALSE])
  weighted_sum <- as.numeric(beta_mat %*% w)
  
  weighted_sum[weighted_sum == 0] <- epsilon
  
  df_weight <- data.frame(
    SNP    = df_coef$SNP,
    weight = 1 / weighted_sum,
    row.names = NULL
  )
  
  list(
    df_coef   = df_coef,
    df_weight = df_weight
  )
}
