# Top-level MIXER function
# ------------------------
#
# Orchestrates:
#   Step 1: SNP ranking by multiple metrics
#   Step 2: Ridge + metric-level weights + SNP-level weights
#   Step 3: Adaptive LASSO selection


MIXER <- function(
    y_train,
    geno_train,
    y_val,
    geno_val,
    top_k = 5000,
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
  # geno_train: training genotype matrix/data.frame (rows: samples, cols: SNPs)
  # y_val: validation outcome vector (0/1)
  # geno_val: validation genotype matrix/data.frame (same cols as geno_train)
  #
  # top_k: number of top SNPs per metric to use for ridge (Step 2)
  # family: glmnet family, default "binomial"
  # nfolds_ridge: CV folds for ridge regression
  # n_threshold: number of thresholds in metric weight calculation
  #
  # min_num: target *minimum* number of selected SNPs for lambda_max (Step 3)
  # max_prop: target *proportion* of SNPs for lambda_min → max_num = max_prop * nrow(df_weight)
  # lambda_init_min, lambda_init_max: initial guesses for lambda range (Step 3)
  # nlambda_adalasso: number of lambda values in final adaptive LASSO grid
  
  # -----------------------------
  # Step 1: SNP ranking by metrics
  # -----------------------------
  message("MIXER - Step 1: SNP ranking by multiple metrics...")
  ranked_list <- rank_snp_metrics(
    phenotype_train = y_train,
    genotype_train  = geno_train,
    phenotype_test  = y_val,
    genotype_test   = geno_val
  )
  
  metric_names <- names(ranked_list)
  if (is.null(metric_names)) {
    stop("rank_snp_metrics must return a *named* list, with one element per metric.")
  }
  
  # Build top-SNP geno lists for train/val
  geno_train_list <- vector("list", length(metric_names))
  geno_val_list   <- vector("list", length(metric_names))
  names(geno_train_list) <- metric_names
  names(geno_val_list)   <- metric_names
  
  for (m in metric_names) {
    ranked_df <- ranked_list[[m]]
    snps_m    <- ranked_df$SNP
    k_m       <- min(top_k, length(snps_m))
    top_snps  <- snps_m[seq_len(k_m)]
    
    common_snps <- intersect(top_snps, colnames(geno_train))
    if (length(common_snps) == 0) {
      stop("No overlapping SNPs between geno_train and ranked SNPs for metric: ", m)
    }
    
    geno_train_list[[m]] <- as.matrix(geno_train[, common_snps, drop = FALSE])
    geno_val_list[[m]]   <- as.matrix(geno_val[,   common_snps, drop = FALSE])
  }
  
  # -----------------------------
  # Step 2a: Ridge per metric
  # -----------------------------
  message("MIXER - Step 2a: Ridge regression for each metric...")
  ridge_coef_list <- lapply(
    geno_train_list,
    function(gm) Ridge_func(
      y      = y_train,
      geno   = gm,
      family = family,
      nfolds = nfolds_ridge
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
    geno_train_list = geno_train_list,
    geno_val_list   = geno_val_list,
    ridge_coef_list = ridge_coef_list,
    n_threshold     = n_threshold
  )
  
  # -----------------------------
  # Step 2c: SNP-level adaptive weights
  # -----------------------------
  message("MIXER - Step 2c: Computing SNP-level adaptive weights...")
  snp_weight_out <- compute_snp_weights(
    ridge_coef_list = ridge_coef_list,
    metric_weights  = metric_weights
  )
  df_coef_all <- snp_weight_out$df_coef
  df_weight   <- snp_weight_out$df_weight
  
  # -----------------------------
  # Step 3: Adaptive LASSO on full genotype
  # -----------------------------
  message("MIXER - Step 3: Running adaptive LASSO...")
  df_selected <- run_adaptive_LASSO(
    y_train        = y_train,
    geno_train     = geno_train,
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
    ranked_snps      = ranked_list,     # Step 1 output
    geno_train_list  = geno_train_list, # top-SNP genotypes per metric (train)
    geno_val_list    = geno_val_list,   # top-SNP genotypes per metric (val)
    ridge_coef_list  = ridge_coef_list, # ridge coefs per metric
    metric_weights   = metric_weights,  # metric-level weights
    snp_coef_all     = df_coef_all,     # union of SNPs with per-metric |beta|
    snp_weights      = df_weight,       # SNP-level adaptive weights
    adaptive_lasso   = df_selected      # final selected SNPs & coefficients
  )
}
