#!/usr/bin/env python3
"""Generate LaTeX rows for the Section 6.2 EOC tables.

The script reads the EOC CSV files produced by
``experiments/section_6_2_convergence/run.py`` and prints the six data rows
for Table 1 (strong imposition) and Table 2 (weak/multiplier imposition).
The theory rows are static in the manuscript and are therefore not generated.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMBINATIONS = (
    ("p1p5", "alpha1p0"),
    ("p1p5", "alpha0p5"),
    ("p2p5", "alpha1p0"),
    ("p2p5", "alpha0p5"),
)
COLUMNS = (
    "eoc_velocity_total",
    "eoc_pressure_pprime",
    "eoc_pressure_l2",
)


def load_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def load_tables(input_dir: Path, enforcement: str):
    tables = []
    for p_code, alpha_code in COMBINATIONS:
        path = input_dir / f"{enforcement}_{p_code}_{alpha_code}_eoc.csv"
        if not path.exists():
            raise SystemExit(f"Missing EOC table: {path}")
        rows = load_rows(path)
        if len(rows) != 6:
            raise SystemExit(f"Expected 6 EOC rows in {path}, found {len(rows)}.")
        tables.append(rows)
    return tables


def format_row(row_index: int, values):
    formatted = [f"{value:.3f}" for value in values]
    line = (
        "\\multicolumn{1}{|c||}{\\cellcolor{lightgray}$"
        f"{row_index + 1}"
        "$} & "
        f"{formatted[0]} & {formatted[1]} & "
        f"\\multicolumn{{1}}{{c||}}{{{formatted[2]}}} & "
        f"{formatted[3]} & {formatted[4]} & "
        f"\\multicolumn{{1}}{{c||}}{{{formatted[5]}}} & "
        f"{formatted[6]} & {formatted[7]} & "
        f"\\multicolumn{{1}}{{c||}}{{{formatted[8]}}} & "
        f"{formatted[9]} & {formatted[10]} & {formatted[11]} \\\\"
    )
    if row_index < 5:
        return line + " \\hline"
    return line + "  \\hline\\hline"


def table_rows(input_dir: Path, enforcement: str):
    tables = load_tables(input_dir, enforcement)
    lines = []
    for index in range(6):
        values = []
        for table in tables:
            row = table[index]
            values.extend(float(row[column]) for column in COLUMNS)
        lines.append(format_row(index, values))
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "published_results" / "section_6_2_convergence",
        help="Directory containing *_eoc.csv files.",
    )
    args = parser.parse_args()

    for title, enforcement in (
        ("Table 1: strong imposition", "strong"),
        ("Table 2: weak/multiplier imposition", "multiplier"),
    ):
        print("% " + title)
        for line in table_rows(args.input, enforcement):
            print("\t\t" + line)
        print()


if __name__ == "__main__":
    main()
