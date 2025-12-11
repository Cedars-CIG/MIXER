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
      feature     = x_train,
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
      feature     = x_train,
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
    feature     = x_train,
    df_weight= df_weight,
    lam_min  = lambda_min,
    lam_max  = lambda_max,
    nlambda  = nlambda,
    family   = family
  )
  
  df_coef
}
