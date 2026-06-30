"""Portable result serialization and EOC summaries."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .metrics import ERROR_KEYS, aggregate_errors, experimental_order, published_scale


def save_result(path, error_series, metadata):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {key: np.asarray(error_series[key], dtype=float) for key in ERROR_KEYS}
    arrays["metadata"] = np.asarray(json.dumps(metadata, sort_keys=True))
    np.savez_compressed(str(path), **arrays)
    return path


def load_result(path):
    with np.load(str(path), allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata"].item()))
        errors = {key: archive[key].astype(float) for key in ERROR_KEYS}
    return errors, metadata


def summarize_result(path):
    errors, metadata = load_result(path)
    aggregates = aggregate_errors(
        errors, metadata["time_step"], metadata["parameters"]["power_law"]
    )
    return {**metadata, **aggregates}


def write_eoc_summary(paths, output_path, minimum_coarse_level=1):
    summaries = sorted(
        (summarize_result(path) for path in paths), key=lambda item: item["level"]
    )
    rows = []
    for coarse, fine in zip(summaries[:-1], summaries[1:]):
        if coarse["level"] < minimum_coarse_level:
            continue
        scale_coarse = published_scale(coarse["mesh_size"], coarse["time_step"])
        scale_fine = published_scale(fine["mesh_size"], fine["time_step"])
        row = {
            "coarse_level": coarse["level"],
            "fine_level": fine["level"],
        }
        for key in ("velocity_total", "pressure_pprime", "pressure_l2"):
            row["eoc_{}".format(key)] = experimental_order(
                coarse[key], fine[key], scale_coarse, scale_fine
            )
        rows.append(row)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "coarse_level",
        "fine_level",
        "eoc_velocity_total",
        "eoc_pressure_pprime",
        "eoc_pressure_l2",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows
