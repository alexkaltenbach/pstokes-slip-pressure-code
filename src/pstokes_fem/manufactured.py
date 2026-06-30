"""Manufactured solution with code generated once by SymPy."""

from __future__ import annotations

from functools import reduce
from math import hypot
from operator import add

import sympy
from fenics import Expression
from scipy import integrate

from .config import time_average_rule


def _ccode(expression) -> str:
    return sympy.printing.ccode(expression).replace("log", "std::log")


def _safe_at_origin(code: str) -> str:
    return (
        "((x[0]*x[0] + x[1]*x[1]) > 1.0e-28 "
        "? ({}) : 0.0)".format(code)
    )


class ManufacturedSolution:
    """Compile all exact fields for fixed spatial and constitutive parameters.

    New ``Expression`` instances only change the scalar parameter ``t``.  The
    generated C++ source is therefore reusable by the legacy FEniCS JIT cache.
    """

    def __init__(self, parameters):
        self.parameters = parameters
        self.degree = parameters.quadrature_degree
        self._build_code()

    def _build_code(self) -> None:
        prm = self.parameters
        time = sympy.symbols("t")
        x0, x1 = sympy.symbols("x[0] x[1]")
        radius = sympy.sqrt(x0**2 + x1**2)

        p_prime = prm.power_law / (prm.power_law - 1.0)
        velocity_exponent = 2.0 * (prm.alpha - 1.0) / prm.power_law + 0.01
        pressure_exponent = prm.alpha - 2.0 / p_prime + 0.01

        velocity = sympy.Matrix(
            [
                time
                * prm.velocity_scale
                * radius**velocity_exponent
                * x1,
                -time
                * prm.velocity_scale
                * radius**velocity_exponent
                * x0,
            ]
        )
        gradient = velocity.jacobian([x0, x1])
        strain = sympy.Rational(1, 2) * (gradient + gradient.T)
        strain_norm = sympy.sqrt(sum(value**2 for value in strain))
        extra_stress = (
            prm.viscosity
            * (prm.delta + strain_norm) ** (prm.power_law - 2.0)
            * strain
        )

        pressure_mean = integrate.dblquad(
            lambda y, x: (
                0.0
                if hypot(x, y) == 0.0 and pressure_exponent < 0.0
                else hypot(x, y) ** pressure_exponent
            ),
            0.0,
            1.0,
            lambda _x: 0.0,
            lambda _x: 1.0,
        )[0]
        pressure = (
            prm.pressure_scale
            * time
            * (radius**pressure_exponent - pressure_mean)
        )

        divergence_stress = sympy.Matrix(
            [
                sympy.diff(extra_stress[0, 0], x0)
                + sympy.diff(extra_stress[0, 1], x1),
                sympy.diff(extra_stress[1, 0], x0)
                + sympy.diff(extra_stress[1, 1], x1),
            ]
        )
        pressure_gradient = sympy.Matrix(
            [sympy.diff(pressure, x0), sympy.diff(pressure, x1)]
        )
        force = sympy.diff(velocity, time) - divergence_stress + pressure_gradient

        self._velocity_code = tuple(_safe_at_origin(_ccode(value)) for value in velocity)
        self._pressure_code = _safe_at_origin(_ccode(pressure))
        self._strain_code = tuple(
            tuple(_safe_at_origin(_ccode(strain[i, j])) for j in range(2))
            for i in range(2)
        )
        self._stress_code = tuple(
            tuple(_safe_at_origin(_ccode(extra_stress[i, j])) for j in range(2))
            for i in range(2)
        )
        self._force_code = tuple(_safe_at_origin(_ccode(value)) for value in force)

    def _expression(self, code, time):
        return Expression(code, t=float(time), degree=self.degree)

    def velocity(self, time):
        return self._expression(self._velocity_code, time)

    def velocity_component(self, component, time):
        if component not in (0, 1):
            raise ValueError("Velocity component must be zero or one.")
        return self._expression(self._velocity_code[component], time)

    def pressure(self, time):
        return self._expression(self._pressure_code, time)

    def strain(self, time):
        return self._expression(self._strain_code, time)

    def stress(self, time):
        return self._expression(self._stress_code, time)

    def force(self, time):
        return self._expression(self._force_code, time)

    def averaged_force(self, t_left, t_right):
        return average_expression(
            self.force, self.parameters.load_rule, t_left, t_right
        )

    def averaged_stress(self, t_left, t_right):
        return average_expression(
            self.stress, self.parameters.load_rule, t_left, t_right
        )


def average_expression(factory, rule, t_left, t_right):
    """Construct a weighted UFL sum of expressions at temporal quadrature nodes."""
    terms = [
        weight * factory(time)
        for weight, time in time_average_rule(rule, t_left, t_right)
    ]
    return reduce(add, terms)

