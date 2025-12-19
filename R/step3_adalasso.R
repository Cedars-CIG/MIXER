# Step 3: Adaptive LASSO selection
# --------------------------------


# Helper 1: given a lambda, how many SNPs are selected?
adaptive_LASSO_selection <- function(y, feature, df_weight, lambda,
                                     family = "binomial") {
  # y: response vector
  # feature: feature matrix/data.frame (rows: samples, cols: SNPs)
  # df_weight: data.frame with columns SNP, weight
  # lambda: scalar penalty value
  # family: glmnet family
  
  x <- as.matrix(feature)
  
  # Ensure SNP order matches df_weight
  common_snps <- intersect(colnames(x), df_weight$SNP)
  if (length(common_snps) == 0) {
    stop("No overlapping SNPs between feature and df_weight in adaptive_LASSO_selection().")
  }
  x <- x[, common_snps, drop = FALSE]
  w <- df_weight$weight[match(common_snps, df_weight$SNP)]
  
  fit <- glmnet::glmnet(
    x      = x,
    y      = y,
    family = family,
    lambda = lambda,
    intercept = TRUE,
    penalty.factor = w
  )
  
  coef_vec <- as.vector(stats::coef(fit))
  # exclude intercept
  support <- which(coef_vec[-1] != 0)
  
  length(support)
}


# Helper 2: full adaptive LASSO fit over a lambda grid, returns selected betas
adaptive_LASSO_fit <- function(y, feature, df_weight,
                               lam_min, lam_max,
                               nlambda = 100,
                               family = "binomial",
                               nfolds = 20) {
  # y: response vector
  # feature: feature matrix/data.frame (rows: samples, cols: SNPs)
  # df_weight: data.frame with columns SNP, weight
  # lam_min, lam_max: lambda range (small = weak penalty, large = strong penalty)
  # nlambda: number of lambda values in the grid
  # family: glmnet family
  # nfolds: CV folds
  
  x <- as.matrix(feature)
  
  # Align SNPs with df_weight
  common_snps <- intersect(colnames(x), df_weight$SNP)
  if (length(common_snps) == 0) {
    stop("No overlapping SNPs between feature and df_weight in adaptive_LASSO_fit().")
  }
  x <- x[, common_snps, drop = FALSE]
  w <- df_weight$weight[match(common_snps, df_weight$SNP)]
  
  # Lambda sequence (log-spaced)
  lam_seq <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))
  
  cvfit <- glmnet::cv.glmnet(
    x      = x,
    y      = y,
    lambda = lam_seq,
    family = family,
    type.measure = "deviance",
    intercept    = TRUE,
    nfolds       = nfolds,
    penalty.factor = w
  )
  
  fit <- glmnet::glmnet(
    x      = x,
    y      = y,
    family = family,
    lambda = cvfit$lambda.min,
    intercept = TRUE,
    penalty.factor = w
  )
  
  message("Step 3 - adaptive LASSO done")
  
  coefs <- stats::coef(fit, s = cvfit$lambda.min)
  coefs <- as.matrix(coefs)
  
  df_coef <- data.frame(
    SNP  = rownames(coefs),
    beta = as.numeric(coefs),
    row.names = NULL
  )
  
  # remove intercept and keep non-zero
  df_coef <- df_coef[df_coef$SNP != "(Intercept)", , drop = FALSE]
  df_coef <- df_coef[df_coef$beta != 0, , drop = FALSE]
  
  df_coef
}


# Main Step-3 wrapper: tune lambda range + run adaptive LASSO

#' Run adaptive LASSO with data-driven lambda range tuning
#'
#' This function tunes \code{lambda_min} and \code{lambda_max} to target a desired
#' range of selected features, then performs cross-validated adaptive LASSO using
#' \code{glmnet} with \code{penalty.factor} determined by \code{df_weight}.
#'
#' @param y_train Training outcome vector.
#' @param feature_train Training feature matrix/data.frame (rows = samples, cols = features).
#' @param df_weight Data.frame with columns \code{SNP} and \code{weight}.
#' @param min_num Target minimum number of selected features when tuning \code{lambda_max}.
#' @param max_prop Target maximum proportion of selected features when tuning \code{lambda_min};
#' \code{max_num = max_prop * nrow(df_weight)}.
#' @param lambda_init_min Initial guess for \code{lambda_min}.
#' @param lambda_init_max Initial guess for \code{lambda_max}.
#' @param nlambda Number of lambdas in the final cross-validated grid.
#' @param max_iter Maximum number of iterations when searching \code{lambda_min} and \code{lambda_max}.
#' @param family glmnet family (default: \code{"binomial"}).
#'
#' @return A data.frame with columns \code{SNP} and \code{beta} for nonzero adaptive LASSO coefficients.
#'
#' @export
run_adaptive_LASSO <- function(
    y_train,
    feature_train,
    df_weight,
    min_num = 10,
    max_prop = 0.05,
    lambda_init_min = 150,
    lambda_init_max = 200,
    nlambda = 100,
    max_iter = 50,
    family = "binomial"
) {
  # y_train: training outcome vector
  # feature_train: training feature matrix/data.frame (rows: samples, cols: SNPs)
  # df_weight: data.frame with columns SNP, weight (adaptive LASSO penalty factors)
  # min_num: target *minimum* number of selected SNPs for lambda_max (upper penalty)
  # max_prop: target *proportion* of SNPs for lambda_min → max_num = max_prop * nrow(df_weight)
  # lambda_init_min, lambda_init_max: initial guesses for lambda range
  # nlambda: number of lambdas for final CV
  # max_iter: max iterations when searching lambda_min and lambda_max
  # family: glmnet family
  
  message("Step 3 - Setting up lambda sequence...")
  
  min_num_target <- min_num
  max_num_target <- max_prop * nrow(df_weight)
  
  lambda_max <- lambda_init_max
  lambda_min <- lambda_init_min
  
  x_train <- as.matrix(feature_train)
  
  common_snps <- intersect(colnames(x_train), df_weight$SNP)
  if (length(common_snps) == 0) {
    stop("No overlapping SNPs between feature_train and df_weight in run_adaptive_LASSO().")
  }
  x_train <- x_train[, common_snps, drop = FALSE]
  
  # Adjust lambda_max: want <= min_num_target SNPs
  num  <- min_num_target + 1
  iter <- 1
  message("Step 3 - Adjusting lambda_max...")
  while ((num > min_num_target) && (iter <= max_iter)) {
    num <- adaptive_LASSO_selection(
      y        = y_train,
      feature  = x_train,
      df_weight= df_weight,
      lambda   = lambda_max,
      family   = family
    )
    
    message("  iter:", iter,
            "  current lambda_max:", lambda_max,
            "  selected SNPs:", num)
    
    if (num > min_num_target) {
      if (num - min_num_target > 20) {
        lambda_max <- 1.2 * lambda_max
      } else {
        lambda_max <- 1.05 * lambda_max
      }
    }
    iter <- iter + 1
  }
  
  # Adjust lambda_min: want >= max_num_target SNPs
  num  <- max_num_target - 1
  iter <- 1
  message("Step 3 - Adjusting lambda_min...")
  while ((num < max_num_target) && (iter <= max_iter)) {
    num <- adaptive_LASSO_selection(
      y        = y_train,
      feature  = x_train,
      df_weight= df_weight,
      lambda   = lambda_min,
      family   = family
    )
    
    message("  iter:", iter,
            "  current lambda_min:", lambda_min,
            "  selected SNPs:", num)
    
    if (num < max_num_target) {
      if (max_num_target - num > 20) {
        lambda_min <- lambda_min / 1.2
      } else {
        lambda_min <- lambda_min / 1.05
      }
    }
    iter <- iter + 1
  }
  
  message("Step 3 - Final lambda_min: ", lambda_min,
          "  lambda_max: ", lambda_max)
  
  # Final adaptive LASSO fit
  df_coef <- adaptive_LASSO_fit(
    y        = y_train,
    feature  = x_train,
    df_weight= df_weight,
    lam_min  = lambda_min,
    lam_max  = lambda_max,
    nlambda  = nlambda,
    family   = family
  )
  
  df_coef
}


#' Evaluate final MIXER model on test data via PRS + threshold moving (best-threshold summary)
#'
#' This function evaluates the final MIXER result on an independent test set using the
#' same procedure as \code{threshold_moving_func()}: compute a PRS on the test set from
#' the final selected coefficients, fit a logistic model \code{y ~ prs}, and run
#' threshold moving to obtain performance metrics across thresholds.
#'
#' In addition to returning the full threshold-moving table, this function also returns
#' a one-row summary at \code{best_threshold}. The \code{best_threshold} is selected as:
#' among thresholds with \code{above_threshold >= prevalence}, choose the threshold whose
#' \code{above_threshold} is closest to the test-set prevalence, where
#' \code{prevalence = mean(y_test)}. If ties occur, the largest threshold is used.
#'
#' @param df_selected Output from \code{run_adaptive_LASSO()}, or a fitted \code{glmnet}/\code{cv.glmnet}
#'   object, or a data.frame/list containing selected features and coefficients.
#' @param y_test Numeric vector (0/1), test outcome.
#' @param feature_test Matrix/data.frame of test features (rows = samples, cols = features).
#' @param threshold_list Optional numeric vector of thresholds passed to \code{threshold_moving_func()}.
#'   If \code{NULL}, thresholds are constructed automatically.
#' @param n_grid Number of thresholds if thresholds are constructed automatically. Default 1000.
#' @param threshold_grid One of \code{c("prob_range", "unit")}.
#'   If \code{"prob_range"}, thresholds are generated from the predicted-probability range.
#'   If \code{"unit"}, thresholds use a fixed grid on 0 to 1.
#' @param s For \code{cv.glmnet} models, which lambda to use (default \code{"lambda.min"}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{prs}: numeric PRS values on the test set
#'   \item \code{glm_prs}: fitted logistic model \code{y ~ prs}
#'   \item \code{threshold_moving}: full data.frame returned by \code{threshold_moving_func()}
#'   \item \code{best_threshold}: numeric best threshold
#'   \item \code{best_row}: one-row data.frame at the best threshold
#'   \item \code{prevalence}: test-set prevalence \code{mean(y_test)}
#' }
#'
#' @export
evaluate_mixer_model <- function(df_selected,
                                 y_test,
                                 feature_test,
                                 threshold_list = NULL,
                                 n_grid = 1000,
                                 threshold_grid = c("prob_range", "unit"),
                                 s = "lambda.min") {
  threshold_grid <- match.arg(threshold_grid)

  feature_test <- as.matrix(feature_test)
  if (nrow(feature_test) != length(y_test)) stop("nrow(feature_test) must equal length(y_test).")
  if (is.null(colnames(feature_test))) stop("feature_test must have column names (feature IDs).")

  # -----------------------------
  # 1) Extract final coefficients
  # -----------------------------
  coef_vec <- NULL

  # Case A: list with $final_model (glmnet or cv.glmnet)
  if (is.list(df_selected) && !is.null(df_selected$final_model) &&
      (inherits(df_selected$final_model, "glmnet") || inherits(df_selected$final_model, "cv.glmnet"))) {

    fm <- df_selected$final_model
    cmat <- if (inherits(fm, "cv.glmnet")) {
      stats::coef(fm, s = s)
    } else {
      s_use <- if (!is.null(fm$lambda) && length(fm$lambda) > 0) fm$lambda[1] else NULL
      stats::coef(fm, s = s_use)
    }

    cmat <- as.matrix(cmat)
    vals <- as.numeric(cmat[, 1])
    names(vals) <- rownames(cmat)
    coef_vec <- vals[setdiff(names(vals), "(Intercept)")]

  # Case B: df_selected itself is glmnet/cv.glmnet
  } else if (inherits(df_selected, "glmnet") || inherits(df_selected, "cv.glmnet")) {

    fm <- df_selected
    cmat <- if (inherits(fm, "cv.glmnet")) {
      stats::coef(fm, s = s)
    } else {
      s_use <- if (!is.null(fm$lambda) && length(fm$lambda) > 0) fm$lambda[1] else NULL
      stats::coef(fm, s = s_use)
    }

    cmat <- as.matrix(cmat)
    vals <- as.numeric(cmat[, 1])
    names(vals) <- rownames(cmat)
    coef_vec <- vals[setdiff(names(vals), "(Intercept)")]

  # Case C: list/data.frame with selected feature + coefficient columns
  } else {
    cand <- NULL
    if (is.list(df_selected) && !is.null(df_selected$df_selected) && is.data.frame(df_selected$df_selected)) {
      cand <- df_selected$df_selected
    } else if (is.list(df_selected) && !is.null(df_selected$selected_df) && is.data.frame(df_selected$selected_df)) {
      cand <- df_selected$selected_df
    } else if (is.data.frame(df_selected)) {
      cand <- df_selected
    }

    if (!is.null(cand)) {
      feat_col <- if ("feature" %in% names(cand)) "feature" else if ("SNP" %in% names(cand)) "SNP" else NA_character_
      coef_col <- if ("coef" %in% names(cand)) "coef" else if ("beta" %in% names(cand)) "beta" else NA_character_
      if (!is.na(feat_col) && !is.na(coef_col)) {
        vals <- cand[[coef_col]]
        names(vals) <- cand[[feat_col]]
        coef_vec <- vals
      }
    }
  }

  if (is.null(coef_vec) || length(coef_vec) == 0) {
    stop("Could not extract final coefficients from df_selected.")
  }

  common_feats <- intersect(names(coef_vec), colnames(feature_test))
  if (length(common_feats) == 0) stop("No overlap between coefficient names and feature_test column names.")
  coef_vec <- coef_vec[common_feats]

  # -----------------------------
  # 2) Compute PRS on test
  # -----------------------------
  prs <- as.numeric(feature_test[, common_feats, drop = FALSE] %*% coef_vec)

  # -----------------------------
  # 3) Fit y ~ prs and run threshold moving
  # -----------------------------
  df_eval <- data.frame(y = y_test, prs = prs)
  glm_prs <- stats::glm(y ~ prs, data = df_eval, family = stats::binomial("logit"))

  # Threshold grid policy (only if user didn't supply threshold_list)
  if (is.null(threshold_list) && threshold_grid == "unit") {
    threshold_list <- seq(0, 1, length.out = n_grid)
  }

  tm <- threshold_moving_func(
    mod = glm_prs,
    df = df_eval,
    threshold_list = threshold_list,
    n_grid = n_grid
  )

  # -----------------------------
  # 4) Best-threshold summary (prevalence-based)
  # -----------------------------
  prevalence <- mean(y_test, na.rm = TRUE)

  df <- tm

  # Keep rows with above_threshold >= prevalence
  df <- df[df$above_threshold >= prevalence, , drop = FALSE]
  if (nrow(df) == 0) {
    warning("No thresholds with above_threshold >= prevalence. best_threshold and best_row set to NULL.")
    best_threshold <- NULL
    best_row <- NULL
  } else {
    df$above_per_distance <- abs(df$above_threshold - prevalence)
    threshold_best <- df$threshold[
      df$above_per_distance == min(df$above_per_distance, na.rm = TRUE)
    ]

    if (length(threshold_best) > 1) {
      threshold_best <- max(threshold_best, na.rm = TRUE)
    }

    best_threshold <- threshold_best
    best_row <- df[df$threshold == threshold_best, , drop = FALSE]
    # Drop helper column if present
    if (!is.null(best_row$above_per_distance)) best_row$above_per_distance <- NULL
  }

  list(
    prs = prs,
    glm_prs = glm_prs,
    threshold_moving = tm,
    best_threshold = best_threshold,
    best_row = best_row,
    prevalence = prevalence
  )
}
