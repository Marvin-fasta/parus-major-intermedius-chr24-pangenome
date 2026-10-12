#!/usr/bin/env python3

import matplotlib.pyplot as plt

# Data: cumulative number of individuals / reference run (13)
n_individuals = list(range(1, 14))

# Structural variants in the local candidate region
structural_variants = [
    164,
    252,
    258,
    266,
    285,
    290,
    285,
    290,
    304,
    305,
    309,
    310,
    431,
]

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

bars = ax.bar(
    n_individuals,
    structural_variants,
    color="#4C78A8",
    edgecolor="black",
    linewidth=0.7,
)

# Axis labels
ax.set_xlabel("Number of individuals", fontsize=14)
ax.set_ylabel("Number of structural variants", fontsize=14)

# Integer ticks for the individual counts
ax.set_xticks(n_individuals)

# Mark run 13 as the reference-inclusive analysis
ax.set_xticklabels(
    [str(i) for i in range(1, 13)] + ["13\n+ reference"],
    fontsize=11,
)

# Thousands separator on y-axis
ax.yaxis.set_major_formatter(
    plt.FuncFormatter(
        lambda value, position: f"{int(value):,}".replace(",", ".")
    )
)

# Print counts above bars
for bar, value in zip(bars, structural_variants):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 5,
        f"{value:,}".replace(",", "."),
        ha="center",
        va="bottom",
        fontsize=9,
        rotation=90,
    )

# Add top margin for labels
ax.set_ylim(0, max(structural_variants) * 1.15)

# Horizontal grid
ax.grid(axis="y", linestyle="--", alpha=0.5)
ax.set_axisbelow(True)

plt.tight_layout()

# Save figures
plt.savefig(
    "local_region_structural_variants_by_number_of_individuals.png",
    dpi=300,
    bbox_inches="tight",
)

plt.savefig(
    "local_region_structural_variants_by_number_of_individuals.pdf",
    bbox_inches="tight",
)

plt.show()
