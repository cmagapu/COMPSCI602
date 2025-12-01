"""
Simplest approach: Monkey-patch DGEB's PairClassificationEvaluator to add ROC-AUC
No custom tasks needed - just patch the original evaluator class
"""

from sklearn.metrics import roc_auc_score
from dgeb.tasks import PairClassificationEvaluator


# Store the original method
_original_compute_metrics = PairClassificationEvaluator._compute_metrics


def _compute_metrics_with_roc_auc(scores, labels, high_score_more_similar):
    """
    Enhanced version of _compute_metrics that adds ROC-AUC.
    Calls the original method and adds our metric.
    """
    # Get all the original metrics
    metrics = _original_compute_metrics(scores, labels, high_score_more_similar)
    
    # Add ROC-AUC
    roc_auc = roc_auc_score(
        labels,
        scores * (1 if high_score_more_similar else -1)
    )
    metrics["roc_auc"] = roc_auc
    
    return metrics


# Monkey-patch the class method
PairClassificationEvaluator._compute_metrics = staticmethod(_compute_metrics_with_roc_auc)

print("✅ PairClassificationEvaluator patched to include ROC-AUC metric")
print("   Now all DGEB tasks using this evaluator will report ROC-AUC!")

