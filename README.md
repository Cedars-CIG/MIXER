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
  3. Predictive Importance Metric (PIM) construction  
  4. Feature-level adaptive weight construction  
  5. Adaptive LASSO selection of the final feature set  

### Data Simulation
- `load(features_example)`
  The dataset `features_example` contains a simulated feature matrix with 5,000 individuals and 16,470 features. The features were derived from a subsample of real-world UK Biobank data and subsequently anonymized by randomly permuting individuals, ensuring that no individual-level information is preserved while retaining realistic marginal distributions, sparsity patterns, and correlation structure across features. This dataset is provided for demonstration and testing purposes only.

- `simulation_data()` 
  Simulates a binary outcome from a high-dimensional feature matrix, performs stratified train/validation/test splitting, and returns MIXER-ready datasets for end-to-end analysis.

### Step 1: Feature Ranking
- `rank_feature_metrics()`  
  Computes per-feature accuracy, balanced accuracy, F1 score, precision, recall, p-value, and ROC-AUC using univariate models, and returns ranked feature lists for each metric.

Below is an example illustrating how different metrics prioritize different features. Each table shows the top-ranked features under a specific metric.

<table>
<tr>
<td>

**Top SNPs by ROC-AUC**

| SNP        | AUC   |
|------------|-------|
| rs229566   | 0.6437 |
| rs5756708  | 0.6402 |
| rs739040   | 0.6370 |
| rs80303    | 0.6358 |
| rs738304   | 0.6347 |
| rs2235334  | 0.6346 |

</td>
<td>

**Top SNPs by −log₁₀(p-value)**

| SNP        | −log₁₀(p) |
|------------|-----------|
| rs1885364  | 14.7645 |
| rs5752613  | 14.2251 |
| rs3747032  | 13.6271 |
| rs3747031  | 13.5874 |
| rs2075455  | 12.9599 |
| rs2034113  | 12.7683 |

</td>
<td>

**Top SNPs by Balanced Accuracy**

| SNP        | Balanced Accuracy |
|------------|-------------------|
| rs738304   | 0.6235 |
| rs715559   | 0.6155 |
| rs715565   | 0.6155 |
| rs229566   | 0.6144 |
| rs2187887  | 0.6100 |
| rs2097465  | 0.6093 |

</td>
</tr>
</table>

### Step 2: Quantify Variable Selection Quality
- `compute_PIM()`  
  Quantifies the predictive importance of each feature subset by aggregating normalized performance metrics across validation thresholds.

- `compute_feature_weights()`  
  Integrates ridge regression coefficients with PIM scores to construct adaptive feature-level weights for downstream selection.

### Step 3: Unify Multiple Selection Criterion
- `run_adaptive_LASSO()`  
  Tunes the adaptive LASSO penalty range and fits the final adaptive LASSO model to select informative features.

### Evaluation Utilities

- `evaluate_mixer_model()`  
  Evaluates the final adaptive LASSO model on an independent test set and reports classification performance metrics including accuracy, balanced accuracy, precision, recall, F1 score, and ROC AUC.


---

## Basic Function Usage

```r
# Example inputs:
#   y_train:         binary outcome (training set)
#   feature_train:   feature matrix for training (samples × features)
#   y_val:           binary outcome (validation set)
#   feature_val:     feature matrix for validation (samples × features)
#   y_test:          binary outcome (test set)
#   feature_test:    feature matrix for test set (samples × features)

library(MIXER)

result <- MIXER(
  y_train           = y_train,
  feature_train    = feature_train,
  y_val             = y_val,
  feature_val      = feature_val,
  y_test            = y_test,
  feature_test     = feature_test,
  top_k             = 5000,     # number of top features used per metric in ridge regression
  n_threshold       = 2000,     # number of thresholds for PIM construction
  min_num           = 10,       # minimum target feature count (lambda_max tuning)
  max_prop          = 0.05,     # maximum proportion of features allowed (lambda_min tuning)
  lambda_init_min   = 150,
  lambda_init_max   = 200,
  nlambda_adalasso  = 100,
  threshold_list    = NULL,     # optional threshold grid for test evaluation
  n_grid_eval       = 1000,     # number of thresholds if grid is constructed automatically
  threshold_grid    = "prob_range"  # or "unit" for a fixed [0,1] grid
)


# Returned object:

# Step 1 outputs
result$ranked_features        # list of ranked features by univariate metrics

# Step 2 outputs
result$ridge_coef_list        # ridge regression coefficients per metric
result$PIM                    # Predictive Importance Metric (metric-level weights)
result$feature_coef_all       # union of features with per-metric ridge coefficients
result$feature_weights        # feature-level adaptive weights

# Step 3 outputs
result$adaptive_lasso         # final adaptive LASSO model / selected features

# Final evaluation
result$test_performance       # predictive performance on the test set

