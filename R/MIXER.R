#' MIXER
#'
#' Multi-metric Integration for eXplanatory and prEdictive Ranking.
#'
#' @param y_train Training outcome vector.
#' @param feature_train Training feature matrix/data.frame (samples x features).
#' @param y_val Validation outcome vector.
#' @param feature_val Validation feature matrix/data.frame (samples x features).
#' @param y_test Test outcome vector.
#' @param feature_test Test feature matrix/data.frame (samples x features).
#' @param top_k Number of top features per metric used for ridge regression.
#' @param family Model family passed to glmnet (default: "binomial").
#' @param nfolds_ridge Number of CV folds for ridge regression.
#' @param n_threshold Number of thresholds for PIM construction.
#' @param min_num Minimum target feature count for lambda_max tuning.
#' @param max_prop Maximum proportion of features allowed for lambda_min tuning.
#' @param lambda_init_min Initial guess for lambda_min.
#' @param lambda_init_max Initial guess for lambda_max.
#' @param nlambda_adalasso Number of lambdas for adaptive LASSO CV grid.
#' @param eval_threshold Classification threshold for computing accuracy/F1/etc on test (default 0.5).
#'
#' @return A list containing ranked features, ridge coefficients, PIM, feature weights,
#' adaptive LASSO selected features, and test-set performance of the final model.
#'
#' @export
MIXER <- function(
    y_train,
    feature_train,
    y_val,
    feature_val,
    y_test,
    feature_test,
    top_k = 5000,
    family = "binomial",
    nfolds_ridge = 20,
    n_threshold = 2000,
    min_num = 10,
    max_prop = 0.05,
    lambda_init_min = 150,
    lambda_init_max = 200,
    nlambda_adalasso = 100,
    eval_threshold = 0.5
) {
  # Coerce to matrices to avoid surprises
  feature_train <- as.matrix(feature_train)
  feature_val   <- as.matrix(feature_val)
  feature_test  <- as.matrix(feature_test)

  if (nrow(feature_train) != length(y_train)) stop("nrow(feature_train) must equal length(y_train).")
  if (nrow(feature_val)   != length(y_val))   stop("nrow(feature_val) must equal length(y_val).")
  if (nrow(feature_test)  != length(y_test))  stop("nrow(feature_test) must equal length(y_test).")

  if (is.null(colnames(feature_train))) stop("feature_train must have column names.")
  if (is.null(colnames(feature_val)))   colnames(feature_val)  <- colnames(feature_train)
  if (is.null(colnames(feature_test)))  colnames(feature_test) <- colnames(feature_train)

  # Ensure validation/test columns align to train by name
  common_val  <- intersect(colnames(feature_train), colnames(feature_val))
  common_test <- intersect(colnames(feature_train), colnames(feature_test))
  if (length(common_val) == 0)  stop("No overlapping feature names between train and validation.")
  if (length(common_test) == 0) stop("No overlapping feature names between train and test.")

  feature_train <- feature_train[, common_val, drop = FALSE]
  feature_val   <- feature_val[,   common_val, drop = FALSE]

  # For test: align to the same columns used by train/val
  common_all <- intersect(colnames(feature_train), colnames(feature_test))
  feature_train <- feature_train[, common_all, drop = FALSE]
  feature_val   <- feature_val[,   common_all, drop = FALSE]
  feature_test  <- feature_test[,  common_all, drop = FALSE]

  # -----------------------------
  # Step 1: feature ranking by metrics
  # -----------------------------
  message("MIXER - Step 1: feature ranking by multiple metrics...")
  ranked_list <- rank_feature_metrics(
    y_train       = y_train,
    feature_train = feature_train,
    y_val         = y_val,
    feature_val   = feature_val
  )

  metric_names <- names(ranked_list)
  if (is.null(metric_names) || length(metric_names) == 0) {
    stop("rank_feature_metrics() must return a *named* list, with one element per metric.")
  }

  # Build top-feature matrices for train/val
  feature_train_list <- vector("list", length(metric_names))
  feature_val_list   <- vector("list", length(metric_names))
  names(feature_train_list) <- metric_names
  names(feature_val_list)   <- metric_names

  for (m in metric_names) {
    ranked_df <- ranked_list[[m]]

    # Support both old ('SNP') and new ('feature') column names
    id_col <- if ("feature" %in% colnames(ranked_df)) "feature" else if ("SNP" %in% colnames(ranked_df)) "SNP" else NA_character_
    if (is.na(id_col)) {
      stop("Each element of ranked_list must contain a column named 'feature' (preferred) or 'SNP'. Missing in metric: ", m)
    }

    feats_m <- ranked_df[[id_col]]
    k_m     <- min(top_k, length(feats_m))
    top_feats <- feats_m[seq_len(k_m)]

    common_feats <- intersect(top_feats, colnames(feature_train))
    if (length(common_feats) == 0) {
      stop("No overlapping features between feature_train and ranked features for metric: ", m)
    }

    feature_train_list[[m]] <- feature_train[, common_feats, drop = FALSE]
    feature_val_list[[m]]   <- feature_val[,   common_feats, drop = FALSE]
  }

  # -----------------------------
  # Step 2a: Ridge per metric
  # -----------------------------
  message("MIXER - Step 2a: Ridge regression for each metric...")
  ridge_coef_list <- lapply(
    feature_train_list,
    function(fm) Ridge_func(
      y       = y_train,
      feature = fm,
      family  = family,
      nfolds  = nfolds_ridge
    )
  )
  names(ridge_coef_list) <- metric_names

  # -----------------------------
  # Step 2b: PIM via PRS performance on validation
  # -----------------------------
  message("MIXER - Step 2b: Computing PIM (Predictive Importance Metric)...")
  pim_df <- compute_PIM(
    y_train            = y_train,
    y_val              = y_val,
    feature_train_list = feature_train_list,
    feature_val_list   = feature_val_list,
    ridge_coef_list    = ridge_coef_list,
    n_threshold        = n_threshold
  )

  # -----------------------------
  # Step 2c: feature-level adaptive weights
  # -----------------------------
  message("MIXER - Step 2c: Computing feature-level adaptive weights...")
  feature_weight_out <- compute_feature_weights(
    ridge_coef_list = ridge_coef_list,
    pim_df          = pim_df
  )
  df_coef_all <- feature_weight_out$df_coef
  df_weight   <- feature_weight_out$df_weight

  # -----------------------------
  # Step 3: Adaptive LASSO on full feature set
  # -----------------------------
  message("MIXER - Step 3: Running adaptive LASSO...")
  df_selected <- run_adaptive_LASSO(
    y_train         = y_train,
    feature_train   = feature_train,
    df_weight       = df_weight,
    min_num         = min_num,
    max_prop        = max_prop,
    lambda_init_min = lambda_init_min,
    lambda_init_max = lambda_init_max,
    nlambda         = nlambda_adalasso,
    family          = family
  )

  # -----------------------------
  # Evaluate final model on TEST data
  # -----------------------------
  # Expect df_selected to contain a fitted glmnet model, or enough information to build one.
  final_model <- NULL
  if (is.list(df_selected) && !is.null(df_selected$final_model)) {
    final_model <- df_selected$final_model
  } else if (inherits(df_selected, "glmnet") || inherits(df_selected, "cv.glmnet")) {
    final_model <- df_selected
  }

  test_performance <- NULL
  if (!is.null(final_model)) {
    # Get predicted probabilities
    prob_test <- tryCatch(
      {
        if (inherits(final_model, "cv.glmnet")) {
          as.numeric(stats::predict(final_model, newx = feature_test, s = "lambda.min", type = "response"))
        } else {
          # For glmnet object without CV, use its internal lambda[which.min] if present; otherwise first lambda
          s_use <- if (!is.null(final_model$lambda) && length(final_model$lambda) > 0) final_model$lambda[1] else NULL
          as.numeric(stats::predict(final_model, newx = feature_test, s = s_use, type = "response"))
        }
      },
      error = function(e) NULL
    )

    if (!is.null(prob_test)) {
      y_pred <- ifelse(prob_test >= eval_threshold, 1, 0)

      # Accuracy
      acc <- mean(y_pred == y_test)

      # AUC (safe)
      auc <- NA_real_
      if (length(unique(y_test)) == 2) {
        pred_obj <- ROCR::prediction(prob_test, y_test)
        auc_obj  <- ROCR::performance(pred_obj, "auc")
        auc <- as.numeric(auc_obj@y.values[[1]])
      }

      # Precision / recall / F1 (safe)
      precision <- recall <- f1 <- NA_real_
      meas <- tryCatch(ROSE::accuracy.meas(y_test, prob_test, threshold = eval_threshold), error = function(e) NULL)
      if (!is.null(meas) && length(meas) >= 5) {
        precision <- meas[[3]]
        recall    <- meas[[4]]
        f1        <- meas[[5]]
      }

      # Balanced accuracy
      df_temp <- data.frame(
        truth     = factor(y_test, levels = c(0, 1)),
        predicted = factor(y_pred, levels = c(0, 1))
      )
      ba_tbl <- yardstick::bal_accuracy(df_temp, truth = truth, estimate = predicted)
      bal_acc <- ba_tbl$.estimate[1]

      test_performance <- list(
        threshold = eval_threshold,
        accuracy = acc,
        balanced_accuracy = bal_acc,
        precision = precision,
        recall = recall,
        f1 = f1,
        roc_auc = auc
      )
    }
  } else {
    warning("Could not locate a fitted final model in adaptive_lasso output; test performance not computed.")
  }

  # -----------------------------
  # Return full result object
  # -----------------------------
  list(
    ranked_features     = ranked_list,          # Step 1 output
    feature_train_list  = feature_train_list,   # top features per metric (train)
    feature_val_list    = feature_val_list,     # top features per metric (val)
    ridge_coef_list     = ridge_coef_list,      # ridge coefs per metric
    PIM                 = pim_df,               # metric-level PIM
    feature_coef_all    = df_coef_all,          # union of features with per-metric |beta|
    feature_weights     = df_weight,            # feature-level adaptive weights
    adaptive_lasso      = df_selected,          # final selected features & coefficients / model
    test_performance    = test_performance      # final model performance on test set
  )
}
