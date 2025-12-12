# Top-level MIXER function
# ------------------------
#
# Orchestrates:
#   Step 1: feature ranking by multiple metrics
#   Step 2: Ridge + metric-level weights + feature-level weights
#   Step 3: Adaptive LASSO selection

#' MIXER
#'
#' Multi-metric Integration for eXplanatory and prEdictive Ranking.
#'
#' @param y_train Training outcome vector.
#' @param feature_train Training feature matrix/data.frame (samples x SNPs).
#' @param y_val Validation outcome vector.
#' @param feature_val Validation feature matrix/data.frame (samples x SNPs).
#' @param top_k Number of top SNPs per metric used for ridge regression.
#' @param family Model family passed to glmnet (default: "binomial").
#' @param nfolds_ridge Number of CV folds for ridge regression.
#' @param n_threshold Number of thresholds for metric weighting.
#' @param min_num Minimum target SNP count for lambda_max tuning.
#' @param max_prop Maximum proportion of SNPs allowed for lambda_min tuning.
#' @param lambda_init_min Initial guess for lambda_min.
#' @param lambda_init_max Initial guess for lambda_max.
#' @param nlambda_adalasso Number of lambdas for adaptive LASSO CV grid.
#'
#' @return A list containing ranked SNPs, ridge coefficients, metric weights,
#' SNP weights, and adaptive LASSO selected SNPs.
#'
#' @export
MIXER <- function(
    y_train,
    feature_train,
    y_val,
    feature_val,
    top_k = 500,
    family = "binomial",
    nfolds_ridge = 20,
    n_threshold = 2000,
    min_num = 10,
    max_prop = 0.05,
    lambda_init_min = 150,
    lambda_init_max = 200,
    nlambda_adalasso = 100
) {
  # y_train: training outcome vector (0/1)
  # feature_train: training feature matrix/data.frame (rows: samples, cols: features)
  # y_val: validation outcome vector (0/1)
  # feature_val: validation feature matrix/data.frame (same cols as feature_train)
  #
  # top_k: number of top features per metric to use for ridge (Step 2)
  # family: glmnet family, default "binomial"
  # nfolds_ridge: CV folds for ridge regression
  # n_threshold: number of thresholds in metric weight calculation
  #
  # min_num: target *minimum* number of selected features for lambda_max (Step 3)
  # max_prop: target *proportion* of features for lambda_min → max_num = max_prop * nrow(df_weight)
  # lambda_init_min, lambda_init_max: initial guesses for lambda range (Step 3)
  # nlambda_adalasso: number of lambda values in final adaptive LASSO grid
  
  # -----------------------------
  # Step 1: feature ranking by metrics
  # -----------------------------
  message("MIXER - Step 1: feature ranking by multiple metrics...")
  ranked_list <- rank_feature_metrics(
    y_train = y_train,
    feature_train  = feature_train,
    y_test  = y_val,
    feature_test   = feature_val
  )
  
  metric_names <- names(ranked_list)
  if (is.null(metric_names)) {
    stop("rank_feature_metrics must return a *named* list, with one element per metric.")
  }
  
  # Build top-feature feature lists for train/val
  feature_train_list <- vector("list", length(metric_names))
  feature_val_list   <- vector("list", length(metric_names))
  names(feature_train_list) <- metric_names
  names(feature_val_list)   <- metric_names
  
  for (m in metric_names) {
    ranked_df <- ranked_list[[m]]
    features_m    <- ranked_df$feature
    k_m       <- min(top_k, length(features_m))
    top_features  <- features_m[seq_len(k_m)]
    
    common_features <- intersect(top_features, colnames(feature_train))
    if (length(common_features) == 0) {
      stop("No overlapping features between feature_train and ranked features for metric: ", m)
    }
    
    feature_train_list[[m]] <- as.matrix(feature_train[, common_features, drop = FALSE])
    feature_val_list[[m]]   <- as.matrix(feature_val[,   common_features, drop = FALSE])
  }
  
  # -----------------------------
  # Step 2a: Ridge per metric
  # -----------------------------
  message("MIXER - Step 2a: Ridge regression for each metric...")
  ridge_coef_list <- lapply(
    feature_train_list,
    function(gm) Ridge_func(
      y       = y_train,
      feature = gm,
      family  = family,
      nfolds  = nfolds_ridge
    )
  )
  names(ridge_coef_list) <- metric_names
  
  # -----------------------------
  # Step 2b: Metric-level weights via PRS performance on validation
  # -----------------------------
  message("MIXER - Step 2b: Computing metric-level weights...")
  metric_weights <- compute_metric_weights(
    y_train         = y_train,
    y_val           = y_val,
    feature_train_list = feature_train_list,
    feature_val_list   = feature_val_list,
    ridge_coef_list = ridge_coef_list,
    n_threshold     = n_threshold
  )
  
  # -----------------------------
  # Step 2c: feature-level adaptive weights
  # -----------------------------
  message("MIXER - Step 2c: Computing feature-level adaptive weights...")
  feature_weight_out <- compute_feature_weights(
    ridge_coef_list = ridge_coef_list,
    metric_weights  = metric_weights
  )
  df_coef_all <- feature_weight_out$df_coef
  df_weight   <- feature_weight_out$df_weight
  
  # -----------------------------
  # Step 3: Adaptive LASSO on full feature
  # -----------------------------
  message("MIXER - Step 3: Running adaptive LASSO...")
  df_selected <- run_adaptive_LASSO(
    y_train        = y_train,
    feature_train     = feature_train,
    df_weight      = df_weight,
    min_num        = min_num,
    max_prop       = max_prop,
    lambda_init_min= lambda_init_min,
    lambda_init_max= lambda_init_max,
    nlambda        = nlambda_adalasso,
    family         = family
  )
  
  # -----------------------------
  # Return full result object
  # -----------------------------
  list(
    ranked_features     = ranked_list,     # Step 1 output
    feature_train_list  = feature_train_list, # top-feature features per metric (train)
    feature_val_list    = feature_val_list,   # top-feature features per metric (val)
    ridge_coef_list  = ridge_coef_list, # ridge coefs per metric
    metric_weights   = metric_weights,  # metric-level weights
    feature_coef_all     = df_coef_all,     # union of features with per-metric |beta|
    feature_weights      = df_weight,       # feature-level adaptive weights
    adaptive_lasso   = df_selected      # final selected features & coefficients
  )
}
