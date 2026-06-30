"""Common UFL operators and the p-Stokes volume residual."""

from fenics import div, dot, dx, grad, inner, sqrt


def symmetric_part(tensor):
    return 0.5 * (tensor + tensor.T)


def tensor_norm(tensor):
    return sqrt(inner(tensor, tensor))


def symmetric_gradient(vector):
    return symmetric_part(grad(vector))


def stress(tensor, power_law, delta, viscosity=1.0):
    return viscosity * (delta + tensor_norm(tensor)) ** (power_law - 2.0) * tensor


def volume_residual(
    velocity,
    pressure,
    pressure_mean_multiplier,
    test_velocity,
    test_pressure,
    test_mean,
    previous_velocity,
    time_step,
    power_law,
    delta,
    viscosity,
    mesh,
):
    """Terms shared by strong and multiplier formulations."""
    strain = symmetric_gradient(velocity)
    test_strain = symmetric_gradient(test_velocity)
    return (
        dot((velocity - previous_velocity) / time_step, test_velocity)
        * dx(domain=mesh)
        + inner(stress(strain, power_law, delta, viscosity), test_strain)
        * dx(domain=mesh)
        - pressure * div(test_velocity) * dx(domain=mesh)
        + div(velocity) * test_pressure * dx(domain=mesh)
        + pressure * test_mean * dx(domain=mesh)
        + pressure_mean_multiplier * test_pressure * dx(domain=mesh)
    )

