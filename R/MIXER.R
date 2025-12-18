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
#' @param eval_threshold Classification threshold for computing test metrics (default 0.5).
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
  # Coerce to matrices
  feature_train <- as.matrix(feature_train)
  feature_val   <- as.matrix(feature_val)
  feature_test  <- as.matrix(feature_test)

  if (nrow(feature_train) != length(y_train)) stop("nrow(feature_train) must equal length(y_train).")
  if (nrow(feature_val)   != length(y_val))   stop("nrow(feature_val) must equal length(y_val).")
  if (nrow(feature_test)  != length(y_test))  stop("nrow(feature_test) must equal length(y_test).")

  if (is.null(colnames(feature_train))) stop("feature_train must have column names.")
  if (is.null(colnames(feature_val)))   colnames(feature_val)  <- colnames(feature_train)
  if (is.null(colnames(feature_test)))  colnames(feature_test) <- colnames(feature_train)

  # Align columns by name
  common <- Reduce(intersect, list(colnames(feature_train), colnames(feature_val), colnames(feature_test)))
  if (length(common) == 0) stop("No overlapping feature names among train/val/test matrices.")
  feature_train <- feature_train[, common, drop = FALSE]
  feature_val   <- feature_val[,   common, drop = FALSE]
  feature_test  <- feature_test[,  common, drop = FALSE]

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
    id_col <- if ("feature" %in% colnames(ranked_df)) "feature" else if ("SNP" %in% colnames(ranked_df)) "SNP" else NA_character_
    if (is.na(id_col)) stop("rank_feature_metrics output must contain 'feature' (preferred) or 'SNP' column. Metric: ", m)

    feats_m <- ranked_df[[id_col]]
    k_m <- min(top_k, length(feats_m))
    top_feats <- feats_m[seq_len(k_m)]

    common_feats <- intersect(top_feats, colnames(feature_train))
    if (length(common_feats) == 0) stop("No overlapping features for metric: ", m)

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
  # Test evaluation (exported helper)
  # -----------------------------
  test_performance <- evaluate_mixer_model(
    df_selected   = df_selected,
    y_test        = y_test,
    feature_test  = feature_test,
    eval_threshold = eval_threshold
  )

  list(
    ranked_features     = ranked_list,
    feature_train_list  = feature_train_list,
    feature_val_list    = feature_val_list,
    ridge_coef_list     = ridge_coef_list,
    PIM                 = pim_df,
    feature_coef_all    = df_coef_all,
    feature_weights     = df_weight,
    adaptive_lasso      = df_selected,
    test_performance    = test_performance
  )
}
