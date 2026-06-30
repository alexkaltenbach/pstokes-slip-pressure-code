"""Command-line entry point for convergence runs."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import (
    ExperimentParameters,
    SolverParameters,
    VALID_ENFORCEMENTS,
    VALID_LOAD_RULES,
)
from .experiment import result_filename, run_level
from .results import write_eoc_summary


def parse_levels(value):
    """Parse `3`, `0:7`, or `0,2,4` into a list of refinement levels."""
    if ":" in value:
        start_text, end_text = value.split(":", 1)
        start, end = int(start_text), int(end_text)
        if end < start:
            raise argparse.ArgumentTypeError("Level range must be increasing.")
        return list(range(start, end + 1))
    if "," in value:
        return [int(item) for item in value.split(",")]
    return [int(value)]


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run the fully discrete p-Stokes convergence experiment."
    )
    parser.add_argument("--enforcement", required=True, choices=VALID_ENFORCEMENTS)
    parser.add_argument("--p", required=True, type=float, dest="power_law")
    parser.add_argument("--alpha", required=True, type=float)
    parser.add_argument("--levels", type=parse_levels, default=parse_levels("0:7"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--load-rule", choices=VALID_LOAD_RULES, default="right")
    parser.add_argument("--pressure-amplitude", type=float)
    parser.add_argument("--absolute-tolerance", type=float, default=1.0e-10)
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--maximum-iterations", type=int, default=50)
    parser.add_argument("--nonlinear-solver", choices=("snes", "newton"), default="snes")
    parser.add_argument("--linear-solver", default="mumps")
    parser.add_argument("--line-search", default="bt")
    parser.add_argument("--report", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv=None):
    arguments = build_parser().parse_args(argv)
    experiment = ExperimentParameters(
        power_law=arguments.power_law,
        alpha=arguments.alpha,
        enforcement=arguments.enforcement,
        load_rule=arguments.load_rule,
        pressure_amplitude=arguments.pressure_amplitude,
    )
    solver = SolverParameters(
        absolute_tolerance=arguments.absolute_tolerance,
        relative_tolerance=arguments.relative_tolerance,
        maximum_iterations=arguments.maximum_iterations,
        linear_solver=arguments.linear_solver,
        nonlinear_solver=arguments.nonlinear_solver,
        line_search=arguments.line_search,
        report=arguments.report,
    )

    arguments.output.mkdir(parents=True, exist_ok=True)
    result_paths = []
    for level in arguments.levels:
        expected = arguments.output / result_filename(experiment, level)
        if expected.exists() and not arguments.overwrite:
            print("Keeping existing {}".format(expected))
            result_paths.append(expected)
            continue
        print(
            "Running enforcement={}, p={}, alpha={}, level={}".format(
                experiment.enforcement,
                experiment.power_law,
                experiment.alpha,
                level,
            )
        )
        result_paths.append(
            run_level(level, experiment, solver, arguments.output)
        )

    if len(result_paths) >= 2:
        summary_name = "{}_p{}_alpha{}_eoc.csv".format(
            experiment.enforcement,
            str(experiment.power_law).replace(".", "p"),
            str(experiment.alpha).replace(".", "p"),
        )
        summary_path = arguments.output / summary_name
        write_eoc_summary(result_paths, summary_path)
        print("Wrote {}".format(summary_path))


if __name__ == "__main__":
    main()
