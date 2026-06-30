"""Small, explicit wrapper around the legacy DOLFIN nonlinear solvers."""

from fenics import (
    Function,
    NewtonSolver,
    NonlinearProblem,
    PETScSNESSolver,
    TrialFunction,
    assemble,
    derivative,
)


class AssembledNonlinearProblem(NonlinearProblem):
    def __init__(self, residual, jacobian, boundary_conditions):
        super(AssembledNonlinearProblem, self).__init__()
        self.residual = residual
        self.jacobian = jacobian
        self.boundary_conditions = list(boundary_conditions)

    def F(self, vector, state):
        assemble(self.residual, tensor=vector)
        for boundary_condition in self.boundary_conditions:
            boundary_condition.apply(vector, state)

    def J(self, matrix, state):
        assemble(self.jacobian, tensor=matrix)
        for boundary_condition in self.boundary_conditions:
            boundary_condition.apply(matrix)


def _make_solver(parameters):
    if parameters.nonlinear_solver == "newton":
        solver = NewtonSolver()
    elif parameters.nonlinear_solver == "snes":
        solver = PETScSNESSolver()
        solver.parameters["method"] = "newtonls"
        solver.parameters["line_search"] = parameters.line_search
    else:
        raise ValueError(
            "Unknown nonlinear solver {!r}.".format(parameters.nonlinear_solver)
        )

    solver.parameters["absolute_tolerance"] = parameters.absolute_tolerance
    solver.parameters["relative_tolerance"] = parameters.relative_tolerance
    solver.parameters["maximum_iterations"] = parameters.maximum_iterations
    solver.parameters["linear_solver"] = parameters.linear_solver
    try:
        solver.parameters["report"] = parameters.report
    except (KeyError, RuntimeError):
        pass
    return solver


def solve_residual(space, residual_builder, boundary_conditions, parameters, guess=None):
    """Build a residual on a coefficient and solve it with a supplied guess."""
    solution = Function(space)
    if guess is not None:
        solution.vector().set_local(guess.vector().get_local())
        solution.vector().apply("insert")

    residual = residual_builder(solution)
    increment = TrialFunction(space)
    jacobian = derivative(residual, solution, increment)
    problem = AssembledNonlinearProblem(residual, jacobian, boundary_conditions)
    solver = _make_solver(parameters)
    solver.solve(problem, solution.vector())
    return solution
