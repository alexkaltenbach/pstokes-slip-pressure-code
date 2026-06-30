"""High-level convergence experiment."""

from __future__ import annotations

from math import sqrt
from pathlib import Path

from fenics import Function, parameters as fenics_parameters, set_log_active

from .errors import ErrorEvaluator
from .manufactured import ManufacturedSolution
from .mesh import unit_square_two_triangles
from .metrics import append_step, empty_error_series
from .problem import solve_time_step
from .results import save_result
from .spaces import build_mixed_space


def configure_fenics():
    fenics_parameters["form_compiler"]["optimize"] = True
    fenics_parameters["form_compiler"]["cpp_optimize"] = True
    fenics_parameters["form_compiler"]["representation"] = "uflacs"
    # Portable default.  Reproducers may opt into -march=native themselves.
    fenics_parameters["form_compiler"]["cpp_optimize_flags"] = "-O3"


def portable_number(value):
    return ("{:g}".format(float(value))).replace("-", "m").replace(".", "p")


def result_filename(parameters, level):
    return "{}_p{}_alpha{}_level{:02d}.npz".format(
        parameters.enforcement,
        portable_number(parameters.power_law),
        portable_number(parameters.alpha),
        int(level),
    )


def run_level(level, parameters, solver_parameters, output_directory):
    configure_fenics()
    set_log_active(bool(solver_parameters.report))

    mesh = unit_square_two_triangles(level)
    space = build_mixed_space(mesh, parameters.enforcement)
    manufactured = ManufacturedSolution(parameters)
    evaluator = ErrorEvaluator(mesh, parameters)
    previous = Function(space)

    number_of_steps = 2 ** (level + 2)
    time_step = parameters.final_time / number_of_steps
    errors = empty_error_series()

    for step in range(1, number_of_steps + 1):
        t_left = (step - 1) * time_step
        t_right = step * time_step
        solution = solve_time_step(
            space,
            previous,
            manufactured,
            t_left,
            t_right,
            parameters,
            solver_parameters,
        )
        components = solution.split(deepcopy=True)
        values = evaluator.evaluate(
            components[0],
            components[1],
            manufactured.velocity(t_right),
            manufactured.pressure(t_right),
            manufactured.strain(t_right),
            manufactured.stress(t_right),
        )
        append_step(errors, values)
        previous.assign(solution)

    metadata = {
        "format_version": 1,
        "level": int(level),
        "mesh_size": sqrt(2.0) * 2.0 ** (-level),
        "time_step": time_step,
        "number_of_steps": number_of_steps,
        "parameters": parameters.as_dict(),
        "solver": solver_parameters.as_dict(),
    }
    path = Path(output_directory) / result_filename(parameters, level)
    return save_result(path, errors, metadata)

