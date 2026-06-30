import pytest


fenics = pytest.importorskip("fenics")

from pstokes_fem.projections import DiscreteLerayProjectors


@pytest.mark.parametrize("enforcement", ("strong", "multiplier"))
def test_discrete_leray_decomposition_and_orthogonality(enforcement):
    mesh = fenics.UnitSquareMesh(1, 1)
    operators = DiscreteLerayProjectors(mesh, enforcement, "default")
    projected, solenoidal, complement = operators.project_basis(0)

    reconstructed = fenics.Function(projected.function_space())
    reconstructed.vector()[:] = solenoidal.vector()[:]
    reconstructed.vector().axpy(1.0, complement.vector())
    reconstructed.vector().axpy(-1.0, projected.vector())

    assert fenics.norm(reconstructed, "L2") <= 1.0e-10
    orthogonality = fenics.assemble(
        fenics.inner(solenoidal, complement) * fenics.dx(domain=mesh)
    )
    assert abs(orthogonality) <= 1.0e-9

