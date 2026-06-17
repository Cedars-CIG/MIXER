#' Simulate outcome data and split into train/validation/test sets
#'
#' This function simulates a binary response \code{y} from a subset of causal
#' features using a logistic model with a calibrated intercept to match a target
#' prevalence. The resulting dataset is split into stratified train, validation,
#' and test sets and is ready for direct use in the MIXER pipeline.
#'
#' If \code{df_feature} is \code{NULL}, the internal example dataset
#' \code{features_example} shipped with the MIXER package is used.
#'
#' @param df_feature A data.frame containing an ID column and feature columns.
#'   If \code{NULL}, uses \code{features_example}.
#' @param effect_size Mean magnitude of causal effects.
#' @param feature_causal_positive Number of causal features with positive effects.
#' @param feature_causal_negative Number of causal features with negative effects.
#' @param target_prev Target prevalence for \code{y = 1}.
#' @param train_prop Proportion of training samples.
#' @param val_prop Proportion of validation samples.
#' @param test_prop Proportion of test samples.
#' @param iid_col Name of the ID column.
#' @param beta_sd Standard deviation around \code{+/-effect_size}.
#' @param seed Random seed for reproducibility.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{y_train}, \code{y_val}, \code{y_test}
#'   \item \code{feature_train}, \code{feature_val}, \code{feature_test}
#'   \item \code{causal_table}
#'   \item \code{intercept}
#'   \item \code{prevalence_observed}
#' }
#'
#' @examples
#' \dontrun{
#' data(features_example)
#' out <- simulation_data()
#' }
#'
#' @export
simulation_data <- function(df_feature = NULL,
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

  # ---- Load example data if needed ----
  if (is.null(df_feature)) {
    utils::data("features_example", package = "MIXER", envir = environment())
    df_feature <- get("features_example", envir = environment(), inherits = FALSE)
  }

  stopifnot(is.data.frame(df_feature))
  if (!iid_col %in% colnames(df_feature)) {
    stop(sprintf("'%s' column not found in df_feature.", iid_col))
  }

  set.seed(seed)

  # ---- Separate IID and feature matrix ----
  feature_names <- setdiff(colnames(df_feature), iid_col)
  X <- as.matrix(df_feature[, feature_names, drop = FALSE])

  n_pos <- feature_causal_positive
  n_neg <- feature_causal_negative
  n_causal <- n_pos + n_neg

  if (n_causal > ncol(X)) {
    stop("Number of causal features exceeds available features.")
  }

  # ---- Simulate causal effects ----
  causal_features <- sample(feature_names, n_causal)
  beta <- c(
    stats::rnorm(n_pos, effect_size, beta_sd),
    stats::rnorm(n_neg, -effect_size, beta_sd)
  )

  causal_table <- data.frame(
    feature = causal_features,
    beta = beta,
    stringsAsFactors = FALSE
  )

  X_causal <- X[, causal_features, drop = FALSE]

  # ---- Simulate outcome ----
  logit_inv <- function(x) 1 / (1 + exp(-x))
  intercept <- stats::uniroot(
    function(b0) mean(logit_inv(drop(X_causal %*% beta) + b0)) - target_prev,
    interval = c(-200, 200)
  )$root

  prob <- logit_inv(drop(X_causal %*% beta) + intercept)
  y <- stats::rbinom(length(prob), 1, prob)

  # ---- Stratified split ----
  idx_case <- which(y == 1)
  idx_ctrl <- which(y == 0)

  split_group <- function(idx) {
    n <- length(idx)
    idx <- sample(idx)
    n_train <- round(train_prop * n)
    n_val <- round(val_prop * n)
    list(
      train = idx[seq_len(n_train)],
      val = idx[(n_train + 1):(n_train + n_val)],
      test = idx[(n_train + n_val + 1):n]
    )
  }

  sc <- split_group(idx_case)
  sn <- split_group(idx_ctrl)

  idx_train <- c(sc$train, sn$train)
  idx_val   <- c(sc$val,   sn$val)
  idx_test  <- c(sc$test,  sn$test)

  # ---- Return MIXER-ready objects (NO IID) ----
  list(
    y_train = y[idx_train],
    y_val = y[idx_val],
    y_test = y[idx_test],
    feature_train = X[idx_train, , drop = FALSE],
    feature_val = X[idx_val, , drop = FALSE],
    feature_test = X[idx_test, , drop = FALSE],
    causal_table = causal_table,
    intercept = intercept,
    prevalence_observed = mean(y)
  )
}
