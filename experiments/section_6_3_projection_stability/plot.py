#!/usr/bin/env python3
"""Plot the Section 6.3 stability data in the original manuscript style.

This script is a light adaptation of the original ``Plot.py`` used for the
paper figure.  The plotting code intentionally keeps the manuscript visual choices
(bmh style, star markers, two-column in-panel legends, large labels).  The only
substantial change is that values are read from ``stability.csv`` instead of
the legacy pickle files.
"""

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
P_CODES = {"15": 1.5, "25": 2.5}
PROJECTOR_LABELS = {
    "P_h": r"$\mathcal{P}_{h_i}$",
    "P_h_perp": r"$\mathcal{P}_{h_i}^\perp$",
}


def read_rows(path, quadrature_degree):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [
            row
            for row in csv.DictReader(stream)
            if int(row["quadrature_degree"]) == quadrature_degree
        ]
    if not rows:
        raise ValueError(
            "No rows with quadrature degree {} in {}".format(
                quadrature_degree,
                path,
            )
        )
    return rows


def values_from_csv(rows, enforcement, p_value, projector):
    """Return the dictionary shape used by the original pickle files."""

    selected = [
        row
        for row in rows
        if row["enforcement"] == enforcement
        and float(row["p"]) == p_value
        and row["projector"] == projector
    ]
    selected.sort(key=lambda row: int(row["mesh_index"]))
    values = {"2": [], "r": [], "r_prime": []}
    h_values = []
    seen_h = set()
    for row in selected:
        h_value = float(row["h"])
        if row["mesh_index"] not in seen_h:
            h_values.append(h_value)
            seen_h.add(row["mesh_index"])
        label = row["exponent_label"]
        if label == "p":
            values["r"].append(float(row["c_stab"]))
        elif label == "p_prime":
            values["r_prime"].append(float(row["c_stab"]))
        else:
            values[label].append(float(row["c_stab"]))
    return h_values, values


def plot_original_column(rows, enforcement, args):
    """Direct adaptation of the original two-row plotting script."""

    fig, ax = plt.subplots(2, 1, figsize=(args.width, args.height))

    for (i, p_code) in enumerate(["15", "25"]):
        p_value = P_CODES[p_code]

        h, values_ph = values_from_csv(rows, enforcement, p_value, "P_h")
        h_perp, values_ph_perp = values_from_csv(
            rows,
            enforcement,
            p_value,
            "P_h_perp",
        )
        if h != h_perp:
            raise ValueError("P_h and P_h_perp h-values do not match.")

        ax[i].semilogx(
            h,
            values_ph["2"],
            c="tab:purple",
            marker="*",
            markersize=args.marker_size,
            ls="solid",
            label=r"$r=2$",
            lw=args.line_width,
        )
        ax[i].semilogx(
            h,
            values_ph["r"],
            c="tab:green",
            marker="*",
            markersize=args.marker_size,
            ls="solid",
            label=r"$r=p$",
            lw=args.line_width,
        )
        ax[i].semilogx(
            h,
            values_ph["r_prime"],
            c="tab:blue",
            marker="*",
            markersize=args.marker_size,
            ls="solid",
            label=r"$r=p^\prime$",
            lw=args.line_width,
        )

        ax[i].semilogx(
            h,
            values_ph_perp["2"],
            c="tab:purple",
            marker="*",
            markersize=args.marker_size,
            ls="dashed",
            label=r"$r=2$",
            lw=args.line_width,
        )
        ax[i].semilogx(
            h,
            values_ph_perp["r"],
            c="tab:green",
            marker="*",
            markersize=args.marker_size,
            ls="dashed",
            label=r"$r=p$",
            lw=args.line_width,
        )
        ax[i].semilogx(
            h,
            values_ph_perp["r_prime"],
            c="tab:blue",
            marker="*",
            markersize=args.marker_size,
            ls="dashed",
            label=r"$r=p^\prime$",
            lw=args.line_width,
        )

        ax[i].set_title(
            r"$\,p={%g}\;$" % p_value,
            bbox={
                "facecolor": "whitesmoke",
                "edgecolor": "black",
                "boxstyle": "round,pad=0.25",
            },
            fontsize=args.title_font_size,
        )

        h_legend, l_legend = ax[i].get_legend_handles_labels()
        handles = [
            plt.plot([], marker="", ls="")[0],
            h_legend[0],
            h_legend[1],
            h_legend[2],
            plt.plot([], marker="", ls="")[0],
            h_legend[3],
            h_legend[4],
            h_legend[5],
        ]
        labels = [
            PROJECTOR_LABELS["P_h"],
            l_legend[0],
            l_legend[1],
            l_legend[2],
            PROJECTOR_LABELS["P_h_perp"],
            l_legend[3],
            l_legend[4],
            l_legend[5],
        ]
        legend = ax[i].legend(
            handles,
            labels,
            ncol=2,
            loc="upper right",
            fontsize=args.legend_font_size,
            handlelength=1.5,
            framealpha=0.9,
            markerscale=1.0,
            labelspacing=0.05,
            handletextpad=0.3,
            borderpad=0.15,
            columnspacing=0.5,
        )
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_linewidth(0.5)

        ax[i].grid(True, which="major", ls="-", color="0.75")
        ax[i].grid(True, which="minor", ls=":", color="0.75")
        if i > 0:
            ax[i].set_xlabel(
                r"$h_i$, $i=1,\ldots,39$",
                fontsize=args.label_font_size,
            )
        ax[i].tick_params(axis="both", labelsize=args.tick_font_size)

    fig.supylabel(
        r"$c^{i}_{\mathrm{stab}}(\mathcal{J}_{h_i})$, "
        r"$\mathcal{J}_{h_i}\in "
        r"\{\mathcal{P}_{h_i},\mathcal{P}_{h_i}^\perp\}$, "
        r"$i=1,\ldots,39$",
        fontsize=args.label_font_size,
    )
    fig.tight_layout()
    return fig


def save_figure(fig, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    png_path = output.with_suffix(".png")
    fig.savefig(png_path, bbox_inches="tight")
    print("Wrote {} and {}".format(output, png_path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT
        / "published_results"
        / "section_6_3_projection_stability"
        / "stability.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "results"
        / "section_6_3_projection_stability"
        / "stability.pdf",
    )
    parser.add_argument(
        "--enforcement",
        choices=("strong", "multiplier"),
        default="strong",
        help="Use 'strong' for the left manuscript panel and 'multiplier' for the right panel.",
    )
    parser.add_argument("--quadrature-degree", type=int, default=6)
    parser.add_argument("--width", type=float, default=8.0)
    parser.add_argument("--height", type=float, default=7.6)
    parser.add_argument("--marker-size", type=float, default=10.0)
    parser.add_argument("--line-width", type=float, default=2.0)
    parser.add_argument("--title-font-size", type=float, default=22.0)
    parser.add_argument("--label-font-size", type=float, default=24.0)
    parser.add_argument("--tick-font-size", type=float, default=24.0)
    parser.add_argument("--legend-font-size", type=float, default=18.0)
    parser.add_argument("--usetex", action="store_true")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    mpl.rcParams["text.usetex"] = args.usetex
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["savefig.dpi"] = 300
    plt.style.use("bmh")

    rows = read_rows(args.input, args.quadrature_degree)
    fig = plot_original_column(rows, args.enforcement, args)
    save_figure(fig, args.output)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
