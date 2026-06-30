"""Per-time-step errors used in the convergence tables."""

from fenics import (
    FiniteElement,
    FunctionSpace,
    TensorElement,
    VectorElement,
    assemble,
    dx,
    inner,
    interpolate,
    sqrt,
)

from .operators import stress, symmetric_gradient, tensor_norm


def _scalar_abs(value):
    return sqrt(value * value)


class ErrorEvaluator:
    def __init__(self, mesh, parameters):
        self.mesh = mesh
        self.parameters = parameters
        cell = mesh.ufl_cell()
        degree = parameters.quadrature_degree
        scheme = parameters.quadrature_rule
        self.tensor_space = FunctionSpace(
            mesh,
            TensorElement(
                "Quadrature", cell, degree=degree, quad_scheme=scheme
            ),
        )
        self.vector_space = FunctionSpace(
            mesh,
            VectorElement(
                "Quadrature", cell, degree=degree, quad_scheme=scheme
            ),
        )
        self.scalar_space = FunctionSpace(
            mesh,
            FiniteElement(
                "Quadrature", cell, degree=degree, quad_scheme=scheme
            ),
        )
        self.measure = dx(domain=mesh, degree=degree, scheme=scheme)

    def evaluate(
        self,
        discrete_velocity,
        discrete_pressure,
        exact_velocity,
        exact_pressure,
        exact_strain,
        exact_stress,
    ):
        prm = self.parameters
        power_law = prm.power_law
        p_prime = power_law / (power_law - 1.0)
        delta = prm.delta

        velocity_exact = interpolate(exact_velocity, self.vector_space)
        pressure_exact = interpolate(exact_pressure, self.scalar_space)
        strain_exact = interpolate(exact_strain, self.tensor_space)
        stress_exact = interpolate(exact_stress, self.tensor_space)

        def natural_map(tensor):
            return (delta + tensor_norm(tensor)) ** (
                (power_law - 2.0) / 2.0
            ) * tensor

        def dual_map(tensor):
            return (delta ** (power_law - 1.0) + tensor_norm(tensor)) ** (
                (p_prime - 2.0) / 2.0
            ) * tensor

        discrete_strain = symmetric_gradient(discrete_velocity)
        discrete_stress = stress(
            discrete_strain, power_law, delta, prm.viscosity
        )
        pressure_error = _scalar_abs(discrete_pressure - pressure_exact)

        values = {
            "natural": assemble(
                inner(
                    natural_map(discrete_strain) - natural_map(strain_exact),
                    natural_map(discrete_strain) - natural_map(strain_exact),
                )
                * self.measure
            ),
            "dual_natural": assemble(
                inner(
                    dual_map(discrete_stress) - dual_map(stress_exact),
                    dual_map(discrete_stress) - dual_map(stress_exact),
                )
                * self.measure
            ),
            "pressure_modular": assemble(
                (
                    (delta + tensor_norm(strain_exact)) ** (power_law - 1.0)
                    + pressure_error
                )
                ** (p_prime - 2.0)
                * pressure_error**2
                * self.measure
            ),
            "pressure_pprime": assemble(
                pressure_error**p_prime * self.measure
            ),
            "pressure_l2": assemble(pressure_error**2 * self.measure),
            "velocity_l2": assemble(
                tensor_norm(discrete_velocity - velocity_exact) ** 2
                * self.measure
            ),
        }
        return values

