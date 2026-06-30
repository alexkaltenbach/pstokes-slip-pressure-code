"""Pure NumPy aggregation of errors and experimental convergence orders."""

from math import log, sqrt

import numpy as np


ERROR_KEYS = (
    "natural",
    "dual_natural",
    "pressure_modular",
    "pressure_pprime",
    "pressure_l2",
    "velocity_l2",
)


def empty_error_series():
    return {key: [] for key in ERROR_KEYS}


def append_step(series, values):
    for key in ERROR_KEYS:
        series[key].append(float(values[key]))


def aggregate_errors(series, time_step, power_law):
    p_prime = power_law / (power_law - 1.0)
    natural = sqrt(time_step * float(np.sum(series["natural"])))
    velocity_linf_l2 = sqrt(float(np.max(series["velocity_l2"])))
    return {
        "natural": natural,
        "dual_natural": sqrt(
            time_step * float(np.sum(series["dual_natural"]))
        ),
        "pressure_modular": sqrt(
            time_step * float(np.sum(series["pressure_modular"]))
        ),
        "pressure_pprime": (
            time_step * float(np.sum(series["pressure_pprime"]))
        )
        ** (1.0 / p_prime),
        "pressure_l2": sqrt(
            time_step * float(np.sum(series["pressure_l2"]))
        ),
        "velocity_linf_l2": velocity_linf_l2,
        "velocity_total": natural + velocity_linf_l2,
    }


def experimental_order(error_coarse, error_fine, scale_coarse, scale_fine):
    if min(error_coarse, error_fine, scale_coarse, scale_fine) <= 0.0:
        raise ValueError("Errors and scales must be positive.")
    return log(error_fine / error_coarse) / log(scale_fine / scale_coarse)


def published_scale(mesh_size, time_step):
    return float(mesh_size) + float(time_step)

