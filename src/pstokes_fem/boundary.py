"""The two explicit implementations of the normal-velocity condition."""

from fenics import (
    CompiledSubDomain,
    DirichletBC,
    FacetNormal,
    dS,
    dot,
    ds,
)


def strong_normal_bcs(space, velocity_component_0, velocity_component_1):
    """Impose inhomogeneous normal data on the axis-aligned unit square."""
    vertical = CompiledSubDomain(
        "on_boundary && (near(x[0], 0.0) || near(x[0], 1.0))"
    )
    horizontal = CompiledSubDomain(
        "on_boundary && (near(x[1], 0.0) || near(x[1], 1.0))"
    )
    return [
        DirichletBC(space.sub(0).sub(0), velocity_component_0, vertical),
        DirichletBC(space.sub(0).sub(1), velocity_component_1, horizontal),
    ]


def multiplier_residual(
    velocity,
    lifted_multiplier,
    test_velocity,
    test_multiplier,
    prescribed_velocity,
    mesh,
):
    """Weak normal constraint for the lifted multiplier.

    The code places the complete known stress vector on the right-hand side.
    Consequently, ``lifted_multiplier · n`` approximates ``-pressure`` rather
    than the physical multiplier from Definition 3.11. Use
    :func:`article_multiplier` to reconstruct the latter.
    """
    normal = FacetNormal(mesh)
    boundary_constraint = (
        dot(velocity, normal) * dot(test_multiplier, normal) * ds(domain=mesh)
        - dot(lifted_multiplier, normal)
        * dot(test_velocity, normal)
        * ds(domain=mesh)
        - dot(prescribed_velocity, normal)
        * dot(test_multiplier, normal)
        * ds(domain=mesh)
    )
    interior_extension = (
        dot(lifted_multiplier("+"), normal("+"))
        * dot(test_multiplier("+"), normal("+"))
        * dS(domain=mesh)
    )
    return boundary_constraint + interior_extension


def article_multiplier(lifted_multiplier, prescribed_stress, mesh):
    """Reconstruct the scalar multiplier used in Definition 3.11.

    If ``ell = lifted_multiplier · n`` denotes the unknown represented by the
    restricted RT element, then

    ``lambda_article = -ell - n · prescribed_stress · n``.

    The supplied stress must be the same endpoint value or temporal average
    that was used in the boundary load of the discrete problem.
    """
    normal = FacetNormal(mesh)
    lifted_normal = dot(lifted_multiplier, normal)
    prescribed_normal_stress = dot(prescribed_stress * normal, normal)
    return -lifted_normal - prescribed_normal_stress
