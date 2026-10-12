#!/usr/bin/env python3

import matplotlib.pyplot as plt

# Daten
n_individuals = list(range(1, 13))

structural_variants = [
    30218,
    58068,
    72542,
    78892,
    83911,
    87869,
    94654,
    98309,
    100793,
    100758,
    102834,
    101506,
]

# Figure erstellen
fig, ax = plt.subplots(figsize=(10, 6))

bars = ax.bar(
    n_individuals,
    structural_variants,
    color="#4C78A8",
    edgecolor="black",
    linewidth=0.7,
)

# Achsenbeschriftungen und Titel
ax.set_xlabel("Number of individuals", fontsize=14)
ax.set_ylabel("Number of structural variants", fontsize=14)


# X-Achse: nur ganze Individuenzahlen
ax.set_xticks(n_individuals)

# Y-Achse mit Tausendertrennzeichen
ax.yaxis.set_major_formatter(
    plt.FuncFormatter(lambda value, position: f"{int(value):,}".replace(",", "."))
)

# Werte oberhalb der Balken anzeigen
for bar, value in zip(bars, structural_variants):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 1200,
        f"{value:,}".replace(",", "."),
        ha="center",
        va="bottom",
        fontsize=8,
        rotation=90,
    )

# Etwas Platz ueber dem hoechsten Balken
ax.set_ylim(0, max(structural_variants) * 1.15)

# Dezente horizontale Hilfslinien
ax.grid(axis="y", linestyle="--", alpha=0.5)
ax.set_axisbelow(True)

plt.tight_layout()

# Grafik speichern
plt.savefig(
    "structural_variants_by_number_of_individuals.png",
    dpi=300,
    bbox_inches="tight",
)

plt.savefig(
    "structural_variants_by_number_of_individuals.pdf",
    bbox_inches="tight",
)

plt.show()
