"""Verify this post's constructed examples; Python standard library only."""
from math import isclose, log
import json


def verify():
    predictions = [1.0, 2.0, 3.0]
    targets = [1.0, 2.0, 3.0]
    paired = [(p - t) ** 2 for p, t in zip(predictions, targets)]
    pairwise = [[(p - t) ** 2 for t in targets] for p in predictions]
    wrong_mse = sum(map(sum, pairwise)) / 9
    right_mse = sum(paired) / 3
    assert isclose(wrong_mse, 4 / 3)
    assert right_mse == 0
    # Independently adjustable predictions minimize the all-pairs objective
    # at the target mean, rather than at their corresponding targets.
    target_mean = sum(targets) / len(targets)
    all_mean_mse = sum((target_mean - t) ** 2 for t in targets) / 3
    assert isclose(all_mean_mse, 2 / 3)
    binary_ce, ten_class_ce = -log(1 / 2), -log(1 / 10)
    for soft_target in ([1.0, 0.0], [0.5, 0.5], [0.2, 0.8]):
        ce = -sum(q * log(1 / 2) for q in soft_target)
        assert isclose(ce, binary_ce)
    assert isclose(binary_ce, 0.6931471805599453)
    assert isclose(ten_class_ce, 2.302585092994046)
    a = (855 + 80) / 1000
    b = (882 + 60) / 1000
    assert isclose(a, 0.935) and isclose(b, 0.942)
    assert isclose((b - a) * 100, 0.7)
    assert isclose((60 / 100 - 80 / 100) * 100, -20)
    return {
        "paired_squared_errors": paired,
        "pairwise_squared_errors": pairwise,
        "paired_mse": right_mse,
        "broadcast_mse": wrong_mse,
        "all_mean_prediction_mse": all_mean_mse,
        "uniform_cross_entropy": {"2_classes": binary_ce, "10_classes": ten_class_ce},
        "aggregate_accuracy": {"model_a": a, "model_b": b},
        "status": "All constructed-example checks passed",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
