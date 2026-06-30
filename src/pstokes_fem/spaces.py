"""Taylor-Hood and boundary-multiplier spaces."""

from fenics import (
    FiniteElement,
    FunctionSpace,
    MixedElement,
    VectorElement,
)


def build_mixed_space(mesh, enforcement: str):
    cell = mesh.ufl_cell()
    velocity = VectorElement("P", cell, 2)
    pressure = FiniteElement("P", cell, 1)
    real = FiniteElement("R", cell, 0)

    if enforcement == "strong":
        element = MixedElement([velocity, pressure, real])
    elif enforcement == "multiplier":
        # RT degree 2 has the P1 normal trace required on each facet.  The
        # restriction discards cell-interior RT degrees of freedom.
        lifted_multiplier = FiniteElement("RT", cell, 2)["facet"]
        element = MixedElement([velocity, pressure, real, lifted_multiplier])
    else:
        raise ValueError("Unknown enforcement {!r}.".format(enforcement))
    return FunctionSpace(mesh, element)


def build_multiplier_trace_space(mesh):
    """Standalone trace space used by the matrix-rank diagnostic."""
    element = FiniteElement("RT", mesh.ufl_cell(), 2)["facet"]
    return FunctionSpace(mesh, element)
