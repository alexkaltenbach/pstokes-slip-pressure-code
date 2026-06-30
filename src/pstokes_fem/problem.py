"""Assembly of one backward-Euler step for either boundary treatment."""

from fenics import FacetNormal, TestFunctions, dot, ds, dx, split

from .boundary import multiplier_residual, strong_normal_bcs
from .nonlinear import solve_residual
from .operators import volume_residual


def solve_time_step(
    space,
    previous_solution,
    manufactured,
    t_left,
    t_right,
    parameters,
    solver_parameters,
):
    """Solve all continuation stages and return the target discrete solution."""
    mesh = space.mesh()
    time_step = t_right - t_left
    previous_velocity = split(previous_solution)[0]
    force = manufactured.averaged_force(t_left, t_right)
    prescribed_stress = manufactured.averaged_stress(t_left, t_right)
    prescribed_velocity = manufactured.velocity(t_right)
    normal = FacetNormal(mesh)

    dx_load = dx(
        domain=mesh,
        degree=parameters.quadrature_degree,
        scheme=parameters.quadrature_rule,
    )
    ds_load = ds(
        domain=mesh,
        degree=parameters.quadrature_degree,
        scheme=parameters.quadrature_rule,
    )

    if parameters.enforcement == "strong":
        boundary_conditions = strong_normal_bcs(
            space,
            manufactured.velocity_component(0, t_right),
            manufactured.velocity_component(1, t_right),
        )
    else:
        boundary_conditions = []

    def residual_builder(stage_power):
        def build(solution):
            fields = split(solution)
            tests = TestFunctions(space)
            velocity, pressure, pressure_mean_multiplier = fields[:3]
            test_velocity, test_pressure, test_mean = tests[:3]

            residual = volume_residual(
                velocity,
                pressure,
                pressure_mean_multiplier,
                test_velocity,
                test_pressure,
                test_mean,
                previous_velocity,
                time_step,
                stage_power,
                parameters.delta,
                parameters.viscosity,
                mesh,
            )
            residual -= dot(force, test_velocity) * dx_load
            residual -= dot(prescribed_stress * normal, test_velocity) * ds_load

            if parameters.enforcement == "multiplier":
                lifted_multiplier = fields[3]
                test_multiplier = tests[3]
                residual += multiplier_residual(
                    velocity,
                    lifted_multiplier,
                    test_velocity,
                    test_multiplier,
                    prescribed_velocity,
                    mesh,
                )
            return residual

        return build

    guess = None
    for stage_power in parameters.continuation_exponents():
        guess = solve_residual(
            space,
            residual_builder(stage_power),
            boundary_conditions,
            solver_parameters,
            guess=guess,
        )
    return guess
