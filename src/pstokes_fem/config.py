"""Configuration objects without a dependency on FEniCS."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import sqrt
from typing import Dict, List, Tuple


VALID_ENFORCEMENTS = ("strong", "multiplier")
VALID_LOAD_RULES = ("right", "midpoint", "gauss2", "gauss3")


def pressure_scale(power_law: float) -> float:
    """Return the pressure amplitude used in the article's experiments."""
    values = {1.5: 1.0e-3, 2.5: 1.0e3}
    for exponent, scale in values.items():
        if abs(float(power_law) - exponent) <= 1.0e-12:
            return scale
    raise ValueError(
        "No published pressure scale is defined for p={!r}; provide one "
        "explicitly.".format(power_law)
    )


def time_average_rule(rule: str, t_left: float, t_right: float) -> List[Tuple[float, float]]:
    """Return ``(weight, time)`` pairs approximating an interval average.

    The weights sum to one.  ``right`` reproduces the supplied implementation;
    the other rules approximate ``dt**-1 * integral(f(t), t_left, t_right)``.
    """
    if not t_right > t_left:
        raise ValueError("The time interval must have positive length.")

    midpoint = 0.5 * (t_left + t_right)
    half_width = 0.5 * (t_right - t_left)
    if rule == "right":
        return [(1.0, t_right)]
    if rule == "midpoint":
        return [(1.0, midpoint)]
    if rule == "gauss2":
        offset = half_width / sqrt(3.0)
        return [(0.5, midpoint - offset), (0.5, midpoint + offset)]
    if rule == "gauss3":
        offset = half_width * sqrt(3.0 / 5.0)
        return [
            (5.0 / 18.0, midpoint - offset),
            (4.0 / 9.0, midpoint),
            (5.0 / 18.0, midpoint + offset),
        ]
    raise ValueError("Unknown load rule {!r}.".format(rule))


@dataclass(frozen=True)
class SolverParameters:
    absolute_tolerance: float = 1.0e-10
    relative_tolerance: float = 1.0e-8
    maximum_iterations: int = 50
    linear_solver: str = "mumps"
    nonlinear_solver: str = "snes"
    line_search: str = "bt"
    report: bool = False

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ExperimentParameters:
    power_law: float
    alpha: float
    enforcement: str
    final_time: float = 0.1
    delta: float = 1.0e-5
    viscosity: float = 1.0
    velocity_scale: float = 1.0
    quadrature_degree: int = 6
    quadrature_rule: str = "default"
    load_rule: str = "right"
    pressure_amplitude: float = None

    def __post_init__(self) -> None:
        if self.enforcement not in VALID_ENFORCEMENTS:
            raise ValueError("Unknown enforcement {!r}.".format(self.enforcement))
        if self.load_rule not in VALID_LOAD_RULES:
            raise ValueError("Unknown load rule {!r}.".format(self.load_rule))
        if self.power_law <= 1.0:
            raise ValueError("The power-law exponent must be greater than one.")
        if self.final_time <= 0.0:
            raise ValueError("The final time must be positive.")

    @property
    def pressure_scale(self) -> float:
        if self.pressure_amplitude is not None:
            return float(self.pressure_amplitude)
        return pressure_scale(self.power_law)

    def continuation_exponents(self) -> Tuple[float, ...]:
        """Continuation path used only to initialize the target solve."""
        target = float(self.power_law)
        if abs(target - 2.0) <= 1.0e-12:
            return (2.0,)

        step = 0.25 if target > 2.0 else -0.25
        path = [2.0]
        current = 2.0 + step
        while (
            (step > 0.0 and current < target - 1.0e-12)
            or (step < 0.0 and current > target + 1.0e-12)
        ):
            path.append(round(current, 12))
            current += step
        path.append(target)

        unique = []
        for exponent in path:
            if not any(abs(existing - exponent) <= 1.0e-12 for existing in unique):
                unique.append(exponent)
        return tuple(unique)

    def as_dict(self) -> Dict[str, object]:
        data = asdict(self)
        data["resolved_pressure_amplitude"] = self.pressure_scale
        data["continuation_exponents"] = self.continuation_exponents()
        return data
