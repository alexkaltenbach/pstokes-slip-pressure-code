"""Meshes used in the convergence experiment."""

from fenics import Mesh, MeshEditor, refine


def unit_square_two_triangles(refinements: int = 0):
    """Create the article's initial mesh and uniformly refine it."""
    if refinements < 0:
        raise ValueError("The refinement level must be non-negative.")

    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, "triangle", 2, 2)
    editor.init_vertices(4)
    editor.init_cells(2)
    editor.add_vertex(0, [0.0, 0.0])
    editor.add_vertex(1, [1.0, 0.0])
    editor.add_vertex(2, [1.0, 1.0])
    editor.add_vertex(3, [0.0, 1.0])
    editor.add_cell(0, [0, 1, 2])
    editor.add_cell(1, [2, 3, 0])
    editor.close()

    for _ in range(refinements):
        mesh = refine(mesh)
    return mesh

