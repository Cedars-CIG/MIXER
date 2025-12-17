#' Simulate a binary response and split data into train/validation/test sets
#'
#' This utility simulates a binary outcome \code{y} from a subset of causal features
#' using a logistic model with a calibrated intercept to match a target prevalence.
#' It then performs a stratified split (by \code{y}) into train/validation/test sets.
#'
#' If \code{df_feature} is \code{NULL}, the function will load the internal example
#' dataset \code{features_example} shipped with the MIXER package.
#'
#' @param df_feature A data.frame containing an ID column and feature columns.
#'   If \code{NULL}, uses \code{features_example}.
#' @param effect_size Mean magnitude of causal effects (positive effects centered at
#'   \code{+effect_size}, negative effects centered at \code{-effect_size}).
#' @param feature_causal_positive Number of causal features with positive effects.
#' @param feature_causal_negative Number of causal features with negative effects.
#' @param target_prev Target prevalence for \code{y = 1}. Default is 0.14.
#' @param train_prop Proportion for training set. Default is 0.8.
#' @param val_prop Proportion for validation set. Default is 0.1.
#' @param test_prop Proportion for test set. Default is 0.1.
#' @param iid_col Name of the ID column. Default is \code{"IID"}.
#' @param beta_sd Standard deviation around \code{+/-effect_size} for causal effects.
#' @param seed Random seed for reproducibility.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{y_train}, \code{y_val}, \code{y_test}
#'   \item \code{feature_train}, \code{feature_val}, \code{feature_test}
#'   \item \code{causal_table} (data.frame with causal features and coefficients)
#'   \item \code{intercept} (calibrated intercept)
#'   \item \code{prevalence_observed} (mean of simulated \code{y})
#' }
#'
#' @examples
#' \dontrun{
#' # Use built-in example feature matrix
#' out <- simulate_and_split_for_mixer(
#'   effect_size = 1,
#'   feature_causal_positive = 50,
#'   feature_causal_negative = 50,
#'   target_prev = 0.14,
#'   seed = 1
#' )
#'
#' # Or provide your own feature data.frame
#' # df_feature <- read.csv("my_features.csv")
#' # out <- simulate_and_split_for_mixer(df_feature)
#' }
#'
#' @export
simulate_and_split_for_mixer <- function(df_feature = NULL,
                                        effect_size = 1,
                                        feature_causal_positive = 50,
                                        feature_causal_negative = 50,
                                        target_prev = 0.14,
                                        train_prop = 0.8,
                                        val_prop = 0.1,
                                        test_prop = 0.1,
                                        iid_col = "IID",
                                        beta_sd = 0.2,
                                        seed = 12345) {
  # ---- Load internal example dataset if needed ----
  if (is.null(df_feature)) {
    # Load into local environment; keep this robust for package usage.
    utils::data("features_example", package = "MIXER", envir = environment())
    if (!exists("features_example", envir = environment(), inherits = FALSE)) {
      stop("Internal dataset 'features_example' could not be loaded. ",
           "Check that data/features_example.rda is present in the installed package.")
    }
    df_feature <- get("features_example", envir = environment(), inherits = FALSE)
  }

  # ---- Basic checks ----
  stopifnot(is.data.frame(df_feature))
  if (!iid_col %in% colnames(df_feature)) {
    stop(sprintf("'%s' column not found in df_feature.", iid_col))
  }
  if (ncol(df_feature) < 2) {
    stop("df_feature must contain an ID column and at least one feature column.")
  }

  props <- c(train_prop, val_prop, test_prop)
  if (any(props <= 0)) stop("train_prop, val_prop, test_prop must be positive.")
  if (abs(sum(props) - 1) > 1e-8) stop("train_prop + val_prop + test_prop must sum to 1.")

  set.seed(seed)

  feature_names <- setdiff(colnames(df_feature), iid_col)
  X <- as.matrix(df_feature[, feature_names, drop = FALSE])
  if (!is.numeric(X)) stop("All feature columns must be numeric.")

  M <- ncol(X)
  n_pos <- as.integer(feature_causal_positive)
  n_neg <- as.integer(feature_causal_negative)
  n_causal <- n_pos + n_neg

  if (n_causal <= 0) stop("Total number of causal features must be > 0.")
  if (n_causal > M) stop("Total number of causal features exceeds number of available features.")

  # ---- Sample causal features and simulate coefficients ----
  causal_features <- sample(feature_names, size = n_causal, replace = FALSE)
  beta_pos <- stats::rnorm(n_pos, mean = effect_size, sd = beta_sd)
  beta_neg <- stats::rnorm(n_neg, mean = -effect_size, sd = beta_sd)
  beta_causal <- c(beta_pos, beta_neg)

  causal_table <- data.frame(
    feature = causal_features,
    beta = beta_causal,
    stringsAsFactors = FALSE
  )

  X_causal <- X[, causal_features, drop = FALSE]

  # ---- Calibrate intercept to match target prevalence ----
  logit_inv <- function(x) 1 / (1 + exp(-x))
  find_intercept <- function(intercept) {
    prob <- logit_inv(drop(X_causal %*% beta_causal) + intercept)
    mean(prob) - target_prev
  }

  intercept <- tryCatch(
    stats::uniroot(function(b0) find_intercept(b0), interval = c(-200, 200))$root,
    error = function(e) {
      stop("Failed to find intercept via uniroot. Try adjusting target_prev or widening the interval.")
    }
  )

  prob <- logit_inv(drop(X_causal %*% beta_causal) + intercept)
  y <- stats::rbinom(n = length(prob), size = 1, prob = prob)
  prevalence_observed <- mean(y)

  # ---- Stratified split by y ----
  split_one_group <- function(indices) {
    m <- length(indices)
    perm <- sample(indices, size = m, replace = FALSE)

    n_train <- as.integer(round(train_prop * m))
    n_val   <- as.integer(round(val_prop * m))
    if (n_train + n_val > m) n_val <- max(0L, m - n_train)
    n_test  <- m - n_train - n_val

    idx_train <- perm[seq_len(n_train)]
    idx_val   <- if (n_val > 0) perm[(n_train + 1):(n_train + n_val)] else integer(0)
    idx_test  <- if (n_test > 0) perm[(n_train + n_val + 1):m] else integer(0)

    list(train = idx_train, val = idx_val, test = idx_test)
  }

  idx_case <- which(y == 1)
  idx_ctrl <- which(y == 0)

  sp_case <- split_one_group(idx_case)
  sp_ctrl <- split_one_group(idx_ctrl)

  idx_train <- sample(c(sp_case$train, sp_ctrl$train))
  idx_val   <- sample(c(sp_case$val,   sp_ctrl$val))
  idx_test  <- sample(c(sp_case$test,  sp_ctrl$test))

  feature_train <- df_feature[idx_train, , drop = FALSE]
  feature_val   <- df_feature[idx_val,   , drop = FALSE]
  feature_test  <- df_feature[idx_test,  , drop = FALSE]

  y_train <- y[idx_train]
  y_val   <- y[idx_val]
  y_test  <- y[idx_test]

  list(
    y_train = y_train,
    y_val = y_val,
    y_test = y_test,
    feature_train = feature_train,
    feature_val = feature_val,
    feature_test = feature_test,
    causal_table = causal_table,
    intercept = intercept,
    prevalence_observed = prevalence_observed
  )
}
