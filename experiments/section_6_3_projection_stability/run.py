#!/usr/bin/env python3
"""Compute the projection-stability indicator reported in Section 6.3."""

import argparse
import csv
from pathlib import Path

from fenics import MPI, UnitSquareMesh, parameters, set_log_level

from pstokes_fem.projections import estimate_stability_constants


ROOT = Path(__file__).resolve().parents[2]
FIELDNAMES = (
    "mesh_index",
    "h",
    "enforcement",
    "p",
    "projector",
    "exponent_label",
    "exponent",
    "quadrature_degree",
    "c_stab",
    "skipped_zero_bases",
)


def parse_integer_selection(value):
    """Parse comma-separated integers and inclusive ranges such as 1:39."""
    selected = []
    for item in value.split(","):
        item = item.strip()
        if ":" in item:
            start, stop = (int(part) for part in item.split(":"))
            selected.extend(range(start, stop + 1))
        else:
            selected.append(int(item))
    values = tuple(sorted(set(selected)))
    if not values or values[0] < 1:
        raise argparse.ArgumentTypeError("mesh indices must be positive")
    return values


def parse_float_selection(value):
    values = tuple(float(item) for item in value.split(","))
    if not values or any(item <= 1.0 for item in values):
        raise argparse.ArgumentTypeError("all p values must be greater than 1")
    return values


def parse_choice_selection(value):
    choices = tuple(item.strip() for item in value.split(",") if item.strip())
    invalid = set(choices) - {"strong", "multiplier"}
    if invalid:
        raise argparse.ArgumentTypeError(
            "unknown enforcement: {}".format(", ".join(sorted(invalid)))
        )
    return choices


def read_rows(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def block_is_complete(rows, mesh_index, enforcement, p_values, degrees):
    actual = {
        (
            float(row["p"]),
            row["projector"],
            row["exponent_label"],
            int(row["quadrature_degree"]),
        )
        for row in rows
        if int(row["mesh_index"]) == mesh_index
        and row["enforcement"] == enforcement
    }
    expected = {
        (p_value, projector, label, degree)
        for p_value in p_values
        for projector in ("P_h", "P_h_perp")
        for label in ("2", "p", "p_prime")
        for degree in degrees
    }
    return expected.issubset(actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mesh-indices",
        type=parse_integer_selection,
        default=parse_integer_selection("1:39"),
    )
    parser.add_argument(
        "--enforcements",
        type=parse_choice_selection,
        default=parse_choice_selection("strong,multiplier"),
    )
    parser.add_argument(
        "--p-values",
        type=parse_float_selection,
        default=parse_float_selection("1.5,2.5"),
    )
    parser.add_argument(
        "--quadrature-degrees",
        type=parse_integer_selection,
        default=parse_integer_selection("6"),
        help="Use 6 for the paper; pass e.g. 6,10,14 for sensitivity checks.",
    )
    parser.add_argument("--linear-solver", default="mumps")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "results"
        / "section_6_3_projection_stability"
        / "stability.csv",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    parameters["form_compiler"]["optimize"] = True
    parameters["form_compiler"]["cpp_optimize"] = True
    set_log_level(30)

    rows = [] if args.overwrite else read_rows(args.output)
    for mesh_index in args.mesh_indices:
        for enforcement in args.enforcements:
            if block_is_complete(
                rows,
                mesh_index,
                enforcement,
                args.p_values,
                args.quadrature_degrees,
            ):
                print(
                    "Keeping mesh {}, enforcement {}".format(
                        mesh_index,
                        enforcement,
                    )
                )
                continue

            print(
                "Running mesh {}, enforcement {}".format(
                    mesh_index,
                    enforcement,
                ),
                flush=True,
            )
            mesh = UnitSquareMesh(mesh_index, mesh_index)
            if MPI.size(mesh.mpi_comm()) != 1:
                raise RuntimeError(
                    "The basis traversal in Section 6.3 must be run in serial."
                )
            new_rows = estimate_stability_constants(
                mesh,
                enforcement,
                p_values=args.p_values,
                quadrature_degrees=args.quadrature_degrees,
                linear_solver=args.linear_solver,
            )
            for row in new_rows:
                row["mesh_index"] = mesh_index
            rows = [
                row
                for row in rows
                if not (
                    int(row["mesh_index"]) == mesh_index
                    and row["enforcement"] == enforcement
                    and float(row["p"]) in args.p_values
                    and int(row["quadrature_degree"])
                    in args.quadrature_degrees
                )
            ]
            rows.extend(new_rows)
            rows.sort(
                key=lambda row: (
                    int(row["mesh_index"]),
                    row["enforcement"],
                    float(row["p"]),
                    row["projector"],
                    row["exponent_label"],
                    int(row["quadrature_degree"]),
                )
            )
            write_rows(args.output, rows)

    print("Wrote {}".format(args.output))


if __name__ == "__main__":
    main()
