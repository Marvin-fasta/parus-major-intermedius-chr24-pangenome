#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot SNP density across a 200 kb local reference region.

The x-axis is converted from local reference coordinates back to
the original Parus major chromosome coordinates.

The original candidate interval is highlighted.
"""

import argparse
import csv

import matplotlib.pyplot as plt


def moving_average(values, window_size):
    """
    Calculate a centered moving average.

    At the beginning and end, the average is calculated using the
    available neighbouring windows.
    """
    smoothed = []
    half_window = window_size // 2

    for index in range(len(values)):
        start = max(0, index - half_window)
        end = min(len(values), index + half_window + 1)

        window_values = values[start:end]
        smoothed.append(sum(window_values) / len(window_values))

    return smoothed


def main():
    parser = argparse.ArgumentParser(
        description="Plot SNP density across a local reference region."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input SNP-density TSV, e.g. snp_density_1kb.tsv",
    )

    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output prefix for PNG and PDF.",
    )

    parser.add_argument(
        "--local-region-start",
        type=int,
        default=1400000,
        help=(
            "First base of the extracted local region in original "
            "reference coordinates. Default: 1400000."
        ),
    )

    parser.add_argument(
        "--candidate-start",
        type=int,
        default=1490771,
        help="Candidate-region start in original reference coordinates.",
    )

    parser.add_argument(
        "--candidate-end",
        type=int,
        default=1525695,
        help="Candidate-region end in original reference coordinates.",
    )

    parser.add_argument(
        "--smoothing-windows",
        type=int,
        default=5,
        help=(
            "Number of 1-kb windows for the moving average. "
            "Default: 5, corresponding to approximately 5 kb."
        ),
    )

    args = parser.parse_args()

    positions = []
    snp_counts = []
    snps_per_kb = []

    with open(args.input, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            local_midpoint = float(row["window_midpoint_bp"])

            # Convert local coordinate to original chromosome coordinate.
            chromosome_midpoint = (
                args.local_region_start - 1 + local_midpoint
            )

            positions.append(chromosome_midpoint)
            snp_counts.append(int(row["snp_count"]))
            snps_per_kb.append(float(row["snps_per_kb"]))

    smoothed_density = moving_average(
        snps_per_kb,
        args.smoothing_windows,
    )

    candidate_values = [
        density
        for position, density in zip(positions, snps_per_kb)
        if args.candidate_start <= position <= args.candidate_end
    ]

    outside_values = [
        density
        for position, density in zip(positions, snps_per_kb)
        if not (args.candidate_start <= position <= args.candidate_end)
    ]

    candidate_mean = (
        sum(candidate_values) / len(candidate_values)
        if candidate_values
        else 0
    )

    outside_mean = (
        sum(outside_values) / len(outside_values)
        if outside_values
        else 0
    )

    fig, ax = plt.subplots(figsize=(14, 6))

    # Highlight the original 34,925 bp candidate region.
    ax.axvspan(
        args.candidate_start,
        args.candidate_end,
        color="#E69F00",
        alpha=0.22,
        label="Candidate region",
        zorder=0,
    )

    # Raw 1-kb values.
    ax.plot(
        positions,
        snps_per_kb,
        color="#56B4E9",
        linewidth=1.0,
        alpha=0.65,
        label="SNP density per 1 kb window",
        zorder=2,
    )

    ax.scatter(
        positions,
        snps_per_kb,
        color="#56B4E9",
        s=12,
        alpha=0.70,
        linewidth=0,
        zorder=3,
    )

    # Smoothed SNP-density line.
    ax.plot(
        positions,
        smoothed_density,
        color="#0072B2",
        linewidth=2.5,
        label=(
            f"Moving average "
            f"({args.smoothing_windows} kb)"
        ),
        zorder=4,
    )

    # Candidate-region boundaries.
    ax.axvline(
        args.candidate_start,
        color="#D55E00",
        linestyle="--",
        linewidth=1.2,
        zorder=1,
    )

    ax.axvline(
        args.candidate_end,
        color="#D55E00",
        linestyle="--",
        linewidth=1.2,
        zorder=1,
    )

    # Region label.
    y_top = max(max(snps_per_kb), max(smoothed_density)) * 1.12

    ax.text(
        (args.candidate_start + args.candidate_end) / 2,
        y_top * 0.98,
        "Candidate region",
        ha="center",
        va="top",
        fontsize=11,
        fontweight="bold",
        color="#8A4B00",
    )

    # Summary box.
    ax.text(
        0.01,
        0.98,
        (
            f"Windows: {len(snps_per_kb)}\n"
            f"Candidate-region mean: {candidate_mean:.2f} SNPs/kb\n"
            f"Outside-region mean: {outside_mean:.2f} SNPs/kb\n"
            f"Difference: {candidate_mean - outside_mean:+.2f} SNPs/kb"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "white",
            "edgecolor": "gray",
            "alpha": 0.92,
        },
    )

    ax.set_xlim(
        args.local_region_start,
        args.local_region_start + 200000,
    )

    ax.set_ylim(0, y_top)

    ax.set_xlabel(
        "Position on chromosome 24 (bp)",
        fontsize=13,
    )

    ax.set_ylabel(
        "SNP density (SNPs per kb)",
        fontsize=13,
    )

    ax.ticklabel_format(
        axis="x",
        style="plain",
        useOffset=False,
    )

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.45,
    )

    ax.set_axisbelow(True)

    ax.legend(
        loc="upper right",
        frameon=True,
        fontsize=10,
    )

    plt.tight_layout()

    plt.savefig(
        f"{args.output_prefix}.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.savefig(
        f"{args.output_prefix}.pdf",
        bbox_inches="tight",
    )

    print(f"Candidate-region mean: {candidate_mean:.2f} SNPs/kb")
    print(f"Outside-region mean: {outside_mean:.2f} SNPs/kb")
    print(
        f"Difference: "
        f"{candidate_mean - outside_mean:+.2f} SNPs/kb"
    )
    print(f"Output: {args.output_prefix}.png")
    print(f"Output: {args.output_prefix}.pdf")


if __name__ == "__main__":
    main()
