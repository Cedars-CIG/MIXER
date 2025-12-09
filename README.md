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


test