#' Example feature matrix for MIXER
#'
#' \code{features_example} is a high-dimensional example feature matrix included
#' to demonstrate the end-to-end workflow of the MIXER package. The dataset
#' contains 5,000 individuals and 16,470 numeric features, along with a unique
#' identifier column (\code{IID}).
#'
#' The feature values were derived from a subsample of real-world UK Biobank data
#' and subsequently anonymized by randomly permuting individuals. This procedure
#' ensures that no individual-level information is preserved while retaining
#' realistic marginal distributions, sparsity patterns, and correlation structure
#' across features. The dataset is provided solely for methodological illustration
#' and testing purposes.
#'
#' @format A data.frame with 5,000 rows and 16,471 columns:
#' \describe{
#'   \item{IID}{Unique individual identifier.}
#'   \item{feature_1, feature_2, \dots}{Numeric feature columns.}
#' }
#'
#' @source
#' UK Biobank (subsampled and anonymized).
"features_example"
