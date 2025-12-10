# MIXER: A Multi-Metric Framework for Robust Feature Selection and Predictive Modeling

High-dimensional biomedical datasets routinely contain sparse signals embedded among vast, correlated features, making variable selection central to building models that generalize. Although significance-based selection is widely used across modalities (e.g., imaging, EHR, multi-omics), statistical significance does not guarantee predictive utility, and vice versa. Yet few methods unify inferential and predictive evidence within a single selection framework. We introduce MIXER (Multi-metric Integration for eXplanatory and prEdictive Ranking), a domain-agnostic approach that integrates multiple selection metrics into one consensus model via adaptive weighting that quantifies each criterion’s contribution.



## Core ideas of MIXER
 <tr>
    <td>
      <img src="/image/github.jpg" alt="overview" width="500">
    </td>
 </tr>

* **Significance ≠ Predictiveness.** What’s “important” for inference is often different from what generalizes for prediction.

* **Multi-metric integration.** Combine p-values with predictive metrics (e.g., AUROC, F1, balanced accuracy) rather than choosing one.

* **One deployable model.** MIXER returns a single consensus model, rather than multiple competing models, for downstream deployment.

## Key features

* Compatible with any selection metric, including p-values, AUROC/AUPRC, F1, balanced accuracy, and more
* Adaptive weighting to quantify the relative utility of each metric
* Merges features selected across metrics into a single consensus model
* Moves beyond p-value–only selection to bridge inference-driven and prediction-driven prioritization


## Installation and Loading

1. To install in R library, use:
     ```ruby
     install.packages("devtools") #if you do not have the devtools package
     devtools::install_github("Cedars-CIG/MIXER")
     ```
2. To load, use:
     ```ruby
     library(MIXER)
     ```

## Tutorial
Please refer to https://cedars-cig.github.io/MIXER/ for a getting-started tutorial. 

## Current Functions

### Top-level Function
- **MIXER()**  
  Runs the full MIXER pipeline, including:
  1. SNP ranking using multiple evaluation metrics  
  2. Ridge regression on top-ranked SNP subsets  
  3. Metric-level weighting via PRS threshold scanning  
  4. SNP-level adaptive weight construction  
  5. Adaptive LASSO selection of final SNP set  

### Step 1: SNP Ranking
- **rank_snp_metrics()**  
  Computes per-SNP accuracy, balanced accuracy, F1, precision, p-value, recall, and ROC-AUC, and returns ranked SNP lists for each metric.

### Step 2: Weight Construction
- **Ridge_func()**  
  Performs ridge regression on a selected SNP subset.

- **threshold_moving_func()**  
  Evaluates PRS predictive performance across probability thresholds.

- **compute_metric_weights()**  
  Builds metric-level weights by combining threshold-based performance profiles.

- **compute_snp_weights()**  
  Integrates ridge coefficients and metric weights to compute adaptive SNP-level weights.

### Step 3: Adaptive LASSO
- **adaptive_LASSO_selection()**  
  Computes the number of selected SNPs at a given λ.

- **adaptive_LASSO_fit()**  
  Runs CV-based adaptive LASSO and returns nonzero coefficients.

- **run_adaptive_LASSO()**  
  Tunes λ range and runs the final adaptive LASSO model.

---

## Basic Function Usage

```r
# Example inputs:
#   y_train:      binary phenotype (training)
#   geno_train:   genotype matrix for training (samples x SNPs)
#   y_val:        binary phenotype (validation)
#   geno_val:     genotype matrix for validation

library(MIXER)

result <- MIXER(
  y_train   = y_train,
  geno_train= geno_train,
  y_val     = y_val,
  geno_val  = geno_val,
  top_k     = 5000,     # number of top SNPs used in ridge per metric
  n_threshold = 2000,   # number of thresholds for metric weighting
  min_num   = 10,       # minimum target SNP count (lambda_max tuning)
  max_prop  = 0.05,     # maximum proportion of SNPs allowed (lambda_min tuning)
  lambda_init_min = 150,
  lambda_init_max = 200,
  nlambda_adalasso = 100
)

# Returned object:
# result$ranked_snps      # list of ranked SNPs by metric
# result$ridge_coef_list  # ridge coefficients per metric
# result$metric_weights   # metric-level weights
# result$snp_weights      # SNP-level adaptive weights
# result$adaptive_lasso   # final selected SNPs and coefficients
