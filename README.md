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

## Current Functions

### Top-level Function
- **MIXER()**  
  Runs the full MIXER pipeline, including:
  1. Feature ranking using multiple evaluation metrics  
  2. Ridge regression on top-ranked feature subsets  
  3. Metric-level weighting via PRS threshold scanning  
  4. Feature-level adaptive weight construction  
  5. Adaptive LASSO selection of final feature set  

### Step 1: Feature Ranking
- **rank_feature_metrics()**  
  Computes per-feature accuracy, balanced accuracy, F1, precision, p-value, recall, and ROC-AUC, and returns ranked feature lists for each metric.

### Step 2: Quantify Variable Selection Quality
- **Ridge_func()**  
  Performs ridge regression on a selected feature subset.

- **threshold_moving_func()**  
  Evaluates PRS predictive performance across probability thresholds.

- **compute_metric_weights()**  
  Builds metric-level weights by combining threshold-based performance profiles.

- **compute_feature_weights()**  
  Integrates ridge coefficients and metric weights to compute adaptive feature-level weights.

### Step 3: Unify Multiple Selection Criterion
- **adaptive_LASSO_selection()**  
  Computes the number of selected features at a given λ.

- **adaptive_LASSO_fit()**  
  Runs CV-based adaptive LASSO and returns nonzero coefficients.

- **run_adaptive_LASSO()**  
  Tunes λ range and runs the final adaptive LASSO model.

---

## Basic Function Usage

```r
# Example inputs:
#   y_train:         binary phenotype (training)
#   feature_train:   feature matrix for training (samples x features)
#   y_val:           binary phenotype (validation)
#   feature_val:     feature matrix for validation

library(MIXER)

result <- MIXER(
  y_train   = y_train,
  feature_train= feature_train,
  y_val     = y_val,
  feature_val  = feature_val,
  top_k     = 5000,     # number of top features used in ridge per metric
  n_threshold = 2000,   # number of thresholds for metric weighting
  min_num   = 10,       # minimum target feature count (lambda_max tuning)
  max_prop  = 0.05,     # maximum proportion of features allowed (lambda_min tuning)
  lambda_init_min = 150,
  lambda_init_max = 200,
  nlambda_adalasso = 100
)

# Returned object:
# result$ranked_features      # list of ranked features by metric
# result$ridge_coef_list      # ridge coefficients per metric
# result$metric_weights       # metric-level weights
# result$feature_weights      # feature-level adaptive weights
# result$adaptive_lasso       # final selected features and coefficients
