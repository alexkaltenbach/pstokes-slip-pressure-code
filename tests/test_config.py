import math

import pytest

from pstokes_fem.config import (
    ExperimentParameters,
    pressure_scale,
    time_average_rule,
)


def test_pressure_scale_matches_experiment():
    assert pressure_scale(1.5) == pytest.approx(1.0e-3)
    assert pressure_scale(2.5) == pytest.approx(1.0e3)
    with pytest.raises(ValueError):
        pressure_scale(2.0)


@pytest.mark.parametrize("rule", ["right", "midpoint", "gauss2", "gauss3"])
def test_time_rule_weights_sum_to_one(rule):
    nodes = time_average_rule(rule, 0.2, 0.7)
    assert sum(weight for weight, _ in nodes) == pytest.approx(1.0)


def test_gauss_rules_integrate_low_order_average():
    left, right = 0.2, 0.7
    exact_linear_average = 0.5 * (left + right)
    for rule in ("midpoint", "gauss2", "gauss3"):
        actual = sum(
            weight * time
            for weight, time in time_average_rule(rule, left, right)
        )
        assert actual == pytest.approx(exact_linear_average)

    exact_quadratic_average = (right**3 - left**3) / (3.0 * (right - left))
    for rule in ("gauss2", "gauss3"):
        actual = sum(
            weight * time**2
            for weight, time in time_average_rule(rule, left, right)
        )
        assert actual == pytest.approx(exact_quadratic_average)


def test_continuation_paths_are_explicit():
    strong = ExperimentParameters(1.5, 1.0, "strong")
    multiplier = ExperimentParameters(2.5, 1.0, "multiplier")
    multiplier_shear_thinning = ExperimentParameters(1.5, 1.0, "multiplier")
    assert strong.continuation_exponents() == (2.0, 1.75, 1.5)
    assert multiplier.continuation_exponents() == (2.0, 2.25, 2.5)
    assert multiplier_shear_thinning.continuation_exponents() == (2.0, 1.75, 1.5)
