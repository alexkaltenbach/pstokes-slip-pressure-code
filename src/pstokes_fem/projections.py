"""Discrete Leray projectors used in the experiment of Section 6.3.

The implementations below are algebraically equivalent to the author-supplied
prototype, but assemble and factorize each left-hand side only once per mesh.
Only right-hand sides change while the nodal basis is traversed.
"""

from fenics import (
    Constant,
    DirichletBC,
    FacetNormal,
    FiniteElement,
    Function,
    FunctionSpace,
    LUSolver,
    Measure,
    MixedElement,
    TrialFunction,
    TrialFunctions,
    TestFunction,
    TestFunctions,
    VectorElement,
    VectorFunctionSpace,
    assemble,
    div,
    dot,
    dx,
    facets,
    grad,
    inner,
    interpolate,
    near,
    MeshFunction,
)


def _vertical_boundary(x, on_boundary):
    return on_boundary and (near(x[0], 0.0) or near(x[0], 1.0))


def _horizontal_boundary(x, on_boundary):
    return on_boundary and (near(x[1], 0.0) or near(x[1], 1.0))


def _normal_bcs(space):
    """Homogeneous normal-trace conditions on the axis-aligned square."""
    return (
        DirichletBC(space.sub(0), Constant(0.0), _vertical_boundary),
        DirichletBC(space.sub(1), Constant(0.0), _horizontal_boundary),
    )


def _mixed_normal_bcs(space, vector_component):
    return (
        DirichletBC(
            space.sub(vector_component).sub(0),
            Constant(0.0),
            _vertical_boundary,
        ),
        DirichletBC(
            space.sub(vector_component).sub(1),
            Constant(0.0),
            _horizontal_boundary,
        ),
    )


def interior_facet_measure(mesh):
    """Return the marked interior-facet measure used for the RT trace lift."""
    topological_dimension = mesh.topology().dim()
    facet_dimension = topological_dimension - 1
    mesh.init(facet_dimension, topological_dimension)
    marker = MeshFunction("size_t", mesh, facet_dimension, 0)
    for facet in facets(mesh):
        if not facet.exterior():
            marker[facet] = 1
    return Measure("dS", domain=mesh, subdomain_data=marker)(1)


class _ReusableLinearProblem:
    """An assembled linear problem with a reusable direct factorization."""

    def __init__(self, bilinear_form, bcs=(), linear_solver="mumps"):
        self.bcs = tuple(bcs)
        self.matrix = assemble(bilinear_form)
        for bc in self.bcs:
            bc.apply(self.matrix)
        # Keeping one solver instance with one unchanged operator lets PETSc
        # retain the direct-solver setup/factorization between RHS solves.
        self.solver = LUSolver(self.matrix, linear_solver)

    def solve(self, linear_form, solution):
        rhs = assemble(linear_form)
        for bc in self.bcs:
            bc.apply(rhs)
        self.solver.solve(solution.vector(), rhs)
        return solution


class _StrongNormalProjection:
    def __init__(self, source_space, linear_solver):
        self.space = source_space
        trial = TrialFunction(self.space)
        self.test = TestFunction(self.space)
        bilinear_form = dot(trial, self.test) * dx
        self.problem = _ReusableLinearProblem(
            bilinear_form,
            _normal_bcs(self.space),
            linear_solver,
        )

    def __call__(self, source):
        target = Function(self.space)
        return self.problem.solve(dot(source, self.test) * dx, target)


class _WeakNormalProjection:
    def __init__(self, source_space, linear_solver):
        self.mesh = source_space.mesh()
        cell = self.mesh.ufl_cell()
        velocity = VectorElement("P", cell, 2)
        multiplier = FiniteElement("RT", cell, 2)["facet"]
        self.space = FunctionSpace(
            self.mesh,
            MixedElement([velocity, multiplier]),
        )
        projected, lifted_multiplier = TrialFunctions(self.space)
        self.velocity_test, multiplier_test = TestFunctions(self.space)
        normal = FacetNormal(self.mesh)
        dS_int = interior_facet_measure(self.mesh)
        bilinear_form = (
            dot(projected, self.velocity_test) * dx
            + dot(projected, normal) * dot(multiplier_test, normal) * Measure(
                "ds", domain=self.mesh
            )
            - dot(lifted_multiplier, normal)
            * dot(self.velocity_test, normal)
            * Measure("ds", domain=self.mesh)
            + dot(lifted_multiplier("+"), normal("+"))
            * dot(multiplier_test("+"), normal("+"))
            * dS_int
        )
        self.problem = _ReusableLinearProblem(
            bilinear_form,
            linear_solver=linear_solver,
        )

    def __call__(self, source):
        mixed_solution = Function(self.space)
        self.problem.solve(
            dot(source, self.velocity_test) * dx,
            mixed_solution,
        )
        projected, _ = mixed_solution.split(deepcopy=True)
        return projected


class _StrongComplementProjection:
    def __init__(self, mesh, linear_solver):
        cell = mesh.ufl_cell()
        pressure = FiniteElement("P", cell, 1)
        velocity = VectorElement("P", cell, 2)
        real = FiniteElement("R", cell, 0)
        self.space = FunctionSpace(
            mesh,
            MixedElement([pressure, velocity, real]),
        )
        potential, discrete_gradient, mean_multiplier = TrialFunctions(
            self.space
        )
        self.pressure_test, velocity_test, mean_test = TestFunctions(self.space)
        bilinear_form = (
            inner(discrete_gradient, grad(self.pressure_test)) * dx
            + potential * mean_test * dx
            + mean_multiplier * self.pressure_test * dx
            + inner(discrete_gradient, velocity_test) * dx
            + potential * div(velocity_test) * dx
        )
        self.problem = _ReusableLinearProblem(
            bilinear_form,
            _mixed_normal_bcs(self.space, 1),
            linear_solver,
        )

    def __call__(self, source):
        mixed_solution = Function(self.space)
        self.problem.solve(
            -div(source) * self.pressure_test * dx,
            mixed_solution,
        )
        _, discrete_gradient, _ = mixed_solution.split(deepcopy=True)
        return discrete_gradient


class _WeakComplementProjection:
    def __init__(self, mesh, linear_solver):
        cell = mesh.ufl_cell()
        pressure = FiniteElement("P", cell, 1)
        velocity = VectorElement("P", cell, 2)
        real = FiniteElement("R", cell, 0)
        multiplier = FiniteElement("RT", cell, 2)["facet"]
        self.space = FunctionSpace(
            mesh,
            MixedElement([pressure, velocity, real, multiplier]),
        )
        (
            potential,
            discrete_gradient,
            mean_multiplier,
            lifted_multiplier,
        ) = TrialFunctions(self.space)
        (
            self.pressure_test,
            velocity_test,
            mean_test,
            multiplier_test,
        ) = TestFunctions(self.space)
        normal = FacetNormal(mesh)
        dS_int = interior_facet_measure(mesh)
        boundary = Measure("ds", domain=mesh)
        bilinear_form = (
            inner(discrete_gradient, grad(self.pressure_test)) * dx
            + dot(discrete_gradient, normal)
            * dot(multiplier_test, normal)
            * boundary
            - dot(lifted_multiplier, normal)
            * dot(velocity_test, normal)
            * boundary
            + dot(lifted_multiplier("+"), normal("+"))
            * dot(multiplier_test("+"), normal("+"))
            * dS_int
            + potential * mean_test * dx
            + mean_multiplier * self.pressure_test * dx
            + inner(discrete_gradient, velocity_test) * dx
            + potential * div(velocity_test) * dx
        )
        self.problem = _ReusableLinearProblem(
            bilinear_form,
            linear_solver=linear_solver,
        )

    def __call__(self, source):
        mixed_solution = Function(self.space)
        self.problem.solve(
            -div(source) * self.pressure_test * dx,
            mixed_solution,
        )
        _, discrete_gradient, _, _ = mixed_solution.split(deepcopy=True)
        return discrete_gradient


class DiscreteLerayProjectors:
    """Project basis functions onto ``V_h`` and apply both Leray parts."""

    def __init__(self, mesh, enforcement, linear_solver="mumps"):
        self.mesh = mesh
        self.source_space = VectorFunctionSpace(mesh, "P", 2)
        if enforcement == "strong":
            self.normal_projection = _StrongNormalProjection(
                self.source_space,
                linear_solver,
            )
            self.complement_projection = _StrongComplementProjection(
                mesh,
                linear_solver,
            )
        elif enforcement == "multiplier":
            self.normal_projection = _WeakNormalProjection(
                self.source_space,
                linear_solver,
            )
            self.complement_projection = _WeakComplementProjection(
                mesh,
                linear_solver,
            )
        else:
            raise ValueError("Unknown enforcement {!r}.".format(enforcement))

    def project_basis(self, dof):
        basis = Function(self.source_space)
        basis.vector()[dof] = 1.0
        projected = self.normal_projection(basis)
        complement = self.complement_projection(projected)
        complement_in_velocity_space = interpolate(
            complement,
            projected.function_space(),
        )
        solenoidal = Function(projected.function_space())
        solenoidal.vector()[:] = projected.vector()[:]
        solenoidal.vector().axpy(-1.0, complement_in_velocity_space.vector())
        return projected, solenoidal, complement_in_velocity_space


def lp_norm(function, exponent, quadrature_degree):
    """Compute the vector-valued L^r norm with explicit quadrature order."""
    mesh = function.function_space().mesh()
    measure = Measure(
        "dx",
        domain=mesh,
        metadata={"quadrature_degree": int(quadrature_degree)},
    )
    magnitude_squared = inner(function, function)
    integral = assemble(magnitude_squared ** (0.5 * exponent) * measure)
    return max(float(integral), 0.0) ** (1.0 / exponent)


def estimate_stability_constants(
    mesh,
    enforcement,
    p_values=(1.5, 2.5),
    quadrature_degrees=(6,),
    linear_solver="mumps",
):
    """Evaluate the basis-function indicator from Section 6.3.

    This computes the maximum over projected nodal basis functions.  It is the
    numerical indicator defined in Section 6.3, not the full operator norm.
    """
    operators = DiscreteLerayProjectors(mesh, enforcement, linear_solver)
    exponent_data = []
    for p_value in p_values:
        exponent_data.extend(
            (
                (float(p_value), "2", 2.0),
                (float(p_value), "p", float(p_value)),
                (
                    float(p_value),
                    "p_prime",
                    float(p_value) / (float(p_value) - 1.0),
                ),
            )
        )

    maxima = {}
    skipped = {}
    for p_value, exponent_label, exponent in exponent_data:
        for degree in quadrature_degrees:
            for projector in ("P_h", "P_h_perp"):
                key = (p_value, exponent_label, int(degree), projector)
                maxima[key] = 0.0
                skipped[key] = 0

    for dof in range(operators.source_space.dim()):
        projected, solenoidal, complement = operators.project_basis(dof)
        denominator_cache = {}
        image_norm_cache = {}
        for p_value, exponent_label, exponent in exponent_data:
            for degree in quadrature_degrees:
                norm_key = (exponent, int(degree))
                if norm_key not in denominator_cache:
                    denominator_cache[norm_key] = lp_norm(
                        projected,
                        exponent,
                        degree,
                    )
                denominator = denominator_cache[norm_key]
                for projector, image in (
                    ("P_h", solenoidal),
                    ("P_h_perp", complement),
                ):
                    key = (p_value, exponent_label, int(degree), projector)
                    if denominator <= 1.0e-14:
                        skipped[key] += 1
                        continue
                    image_key = (projector, exponent, int(degree))
                    if image_key not in image_norm_cache:
                        image_norm_cache[image_key] = lp_norm(
                            image,
                            exponent,
                            degree,
                        )
                    ratio = image_norm_cache[image_key] / denominator
                    maxima[key] = max(maxima[key], ratio)

    rows = []
    for p_value, exponent_label, exponent in exponent_data:
        for degree in quadrature_degrees:
            for projector in ("P_h", "P_h_perp"):
                key = (p_value, exponent_label, int(degree), projector)
                rows.append(
                    {
                        "mesh_index": None,
                        "h": float(mesh.hmax()),
                        "enforcement": enforcement,
                        "p": p_value,
                        "projector": projector,
                        "exponent_label": exponent_label,
                        "exponent": exponent,
                        "quadrature_degree": int(degree),
                        "c_stab": maxima[key],
                        "skipped_zero_bases": skipped[key],
                    }
                )
    return rows
