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



#' Evaluate a fitted final model on test data
#'
#' This helper extracts the fitted final model from an adaptive LASSO output
#' (or accepts a glmnet/cv.glmnet model directly), generates predicted
#' probabilities on the test set, and computes standard classification metrics.
#'
#' @param df_selected Output from \code{run_adaptive_LASSO()} or a fitted
#'   \code{glmnet}/\code{cv.glmnet} object. If a list is provided, it should
#'   contain a fitted model in \code{$final_model}.
#' @param y_test Numeric vector (0/1), test outcome.
#' @param feature_test Matrix/data.frame of test features (rows = samples, cols = features).
#' @param eval_threshold Classification threshold for hard predictions. Default 0.5.
#' @param s For \code{cv.glmnet} models, which lambda to use (default \code{"lambda.min"}).
#'   Ignored for plain \code{glmnet} objects.
#'
#' @return A list with elements \code{threshold}, \code{accuracy},
#'   \code{balanced_accuracy}, \code{precision}, \code{recall}, \code{f1},
#'   \code{roc_auc}. Returns \code{NULL} if predictions cannot be computed.
#'
#' @examples
#' \dontrun{
#' out <- simulation_data(seed = 1)
#' res <- MIXER(out$y_train, out$feature_train, out$y_val, out$feature_val,
#'              out$y_test, out$feature_test)
#' perf <- evaluate_mixer_model(res$adaptive_lasso, out$y_test, out$feature_test)
#' }
#'
#' @export
evaluate_mixer_model <- function(df_selected,
                                 y_test,
                                 feature_test,
                                 eval_threshold = 0.5,
                                 s = "lambda.min") {
  feature_test <- as.matrix(feature_test)
  if (nrow(feature_test) != length(y_test)) stop("nrow(feature_test) must equal length(y_test).")

  # Extract model
  final_model <- NULL
  if (is.list(df_selected) && !is.null(df_selected$final_model)) {
    final_model <- df_selected$final_model
  } else if (inherits(df_selected, "glmnet") || inherits(df_selected, "cv.glmnet")) {
    final_model <- df_selected
  }

  if (is.null(final_model)) {
    warning("Could not locate a fitted final model in adaptive_lasso output; test performance not computed.")
    return(NULL)
  }

  # Predicted probabilities
  prob_test <- tryCatch(
    {
      if (inherits(final_model, "cv.glmnet")) {
        as.numeric(stats::predict(final_model, newx = feature_test, s = s, type = "response"))
      } else {
        # For glmnet object without CV, use first lambda by default
        s_use <- if (!is.null(final_model$lambda) && length(final_model$lambda) > 0) final_model$lambda[1] else NULL
        as.numeric(stats::predict(final_model, newx = feature_test, s = s_use, type = "response"))
      }
    },
    error = function(e) NULL
  )

  if (is.null(prob_test)) return(NULL)

  y_pred <- ifelse(prob_test >= eval_threshold, 1, 0)

  # Accuracy
  acc <- mean(y_pred == y_test)

  # ROC AUC (safe)
  auc <- NA_real_
  if (length(unique(y_test)) == 2) {
    pred_obj <- ROCR::prediction(prob_test, y_test)
    auc_obj  <- ROCR::performance(pred_obj, "auc")
    auc <- as.numeric(auc_obj@y.values[[1]])
  }

  # Precision / recall / F1 (safe)
  precision <- recall <- f1 <- NA_real_
  meas <- tryCatch(
    ROSE::accuracy.meas(y_test, prob_test, threshold = eval_threshold),
    error = function(e) NULL
  )
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

  list(
    threshold = eval_threshold,
    accuracy = acc,
    balanced_accuracy = bal_acc,
    precision = precision,
    recall = recall,
    f1 = f1,
    roc_auc = auc
  )
}