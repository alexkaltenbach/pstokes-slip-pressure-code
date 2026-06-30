"""Merged finite-element implementation for the unsteady p-Stokes system."""

from .config import (
    ExperimentParameters,
    SolverParameters,
    pressure_scale,
    time_average_rule,
)

__all__ = [
    "ExperimentParameters",
    "SolverParameters",
    "pressure_scale",
    "time_average_rule",
]

__version__ = "1.0.0"
