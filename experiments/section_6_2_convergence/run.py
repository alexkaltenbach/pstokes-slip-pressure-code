#!/usr/bin/env python3
"""Run every parameter combination used in Tables 1 and 2."""

import argparse
from pathlib import Path

from pstokes_fem.cli import parse_levels
from pstokes_fem.config import ExperimentParameters, SolverParameters
from pstokes_fem.experiment import result_filename, run_level
from pstokes_fem.results import write_eoc_summary


ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "section_6_2_convergence",
    )
    parser.add_argument("--levels", type=parse_levels, default=parse_levels("0:7"))
    parser.add_argument(
        "--load-rule",
        choices=("right", "midpoint", "gauss2", "gauss3"),
        default="right",
    )
    parser.add_argument("--absolute-tolerance", type=float, default=1.0e-10)
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--maximum-iterations", type=int, default=50)
    parser.add_argument("--linear-solver", default="mumps")
    parser.add_argument("--nonlinear-solver", choices=("snes", "newton"), default="snes")
    parser.add_argument("--line-search", default="bt")
    parser.add_argument("--report", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    solver = SolverParameters(
        absolute_tolerance=args.absolute_tolerance,
        relative_tolerance=args.relative_tolerance,
        maximum_iterations=args.maximum_iterations,
        linear_solver=args.linear_solver,
        nonlinear_solver=args.nonlinear_solver,
        line_search=args.line_search,
        report=args.report,
    )
    for enforcement in ("strong", "multiplier"):
        for power_law in (1.5, 2.5):
            for alpha in (1.0, 0.5):
                prm = ExperimentParameters(
                    power_law=power_law,
                    alpha=alpha,
                    enforcement=enforcement,
                    load_rule=args.load_rule,
                )
                paths = []
                for level in args.levels:
                    path = args.output / result_filename(prm, level)
                    if path.exists() and not args.overwrite:
                        print("Keeping existing {}".format(path))
                    else:
                        print(
                            "Running {}, p={}, alpha={}, level={}".format(
                                enforcement, power_law, alpha, level
                            )
                        )
                        path = run_level(level, prm, solver, args.output)
                    paths.append(path)
                if len(paths) >= 2:
                    summary = args.output / (
                        "{}_p{}_alpha{}_eoc.csv".format(
                            enforcement,
                            str(power_law).replace(".", "p"),
                            str(alpha).replace(".", "p"),
                        )
                    )
                    write_eoc_summary(paths, summary)


if __name__ == "__main__":
    main()
