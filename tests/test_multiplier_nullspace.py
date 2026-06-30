"""Rank diagnostic for the facet-restricted RT multiplier."""

import numpy as np
import pytest

pytest.importorskip("fenics")

from fenics import (  # noqa: E402
    CompiledSubDomain,
    FacetNormal,
    Measure,
    MeshFunction,
    TestFunction,
    TrialFunction,
    Vertex,
    assemble,
    dS,
    dot,
    ds,
    facets,
)

from pstokes_fem.mesh import unit_square_two_triangles  # noqa: E402
from pstokes_fem.spaces import build_multiplier_trace_space  # noqa: E402


def marked_measure(mesh):
    mesh.init(1, 0)
    mesh.init(0, 1)
    boundary = MeshFunction("bool", mesh, mesh.topology().dim() - 1, False)
    CompiledSubDomain("on_boundary").mark(boundary, True)
    marker = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)

    for facet in facets(mesh):
        if facet.exterior():
            continue
        touches_boundary = False
        for vertex_index in facet.entities(0):
            vertex = Vertex(mesh, vertex_index)
            if any(boundary[index] for index in vertex.entities(1)):
                touches_boundary = True
                break
        if not touches_boundary:
            marker[facet] = 1
    return Measure("dS", domain=mesh, subdomain_data=marker)


def mass_matrix(mesh, interior_measure):
    space = build_multiplier_trace_space(mesh)
    trial = TrialFunction(space)
    test = TestFunction(space)
    normal = FacetNormal(mesh)
    form = dot(trial, normal) * dot(test, normal) * ds(domain=mesh)
    form += (
        dot(trial("+"), normal("+"))
        * dot(test("+"), normal("+"))
        * interior_measure
    )
    return assemble(form).array()


def numerical_nullity(matrix):
    values = np.linalg.svd(matrix, compute_uv=False)
    tolerance = max(matrix.shape) * np.finfo(float).eps * values[0]
    return int(np.count_nonzero(values <= tolerance))


@pytest.mark.parametrize("level", [0, 1, 2])
def test_original_measure_equals_plain_dS_and_has_full_rank(level):
    mesh = unit_square_two_triangles(level)
    original = mass_matrix(mesh, marked_measure(mesh))
    plain = mass_matrix(mesh, dS(domain=mesh))
    np.testing.assert_allclose(original, plain, rtol=1.0e-12, atol=1.0e-12)
    assert numerical_nullity(original) == 0

