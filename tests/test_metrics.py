import math

import pytest

from pstokes_fem.metrics import aggregate_errors, experimental_order


def test_velocity_uses_linf_in_time():
    series = {
        "natural": [1.0, 4.0],
        "dual_natural": [1.0, 1.0],
        "pressure_modular": [1.0, 1.0],
        "pressure_pprime": [1.0, 1.0],
        "pressure_l2": [1.0, 1.0],
        "velocity_l2": [1.0, 9.0],
    }
    values = aggregate_errors(series, time_step=0.5, power_law=1.5)
    assert values["natural"] == pytest.approx(math.sqrt(2.5))
    assert values["velocity_linf_l2"] == pytest.approx(3.0)
    assert values["velocity_total"] == pytest.approx(math.sqrt(2.5) + 3.0)


def test_eoc_for_halved_scale():
    assert experimental_order(0.25, 0.125, 1.0, 0.5) == pytest.approx(1.0)

