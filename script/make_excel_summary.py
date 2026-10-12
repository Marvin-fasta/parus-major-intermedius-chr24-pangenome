#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

import pandas as pd


results_dir = Path(sys.argv[1])
samples = sys.argv[2].split(",")

rows = []

for run_number in range(1, len(samples) + 1):
    run_id = f"{run_number:02d}"
    run_dir = results_dir / "runs" / f"run_{run_id}"

    stats_file = run_dir / "variants.stats.txt"
    sv_summary_file = run_dir / "sv_summary.tsv"

    if not stats_file.exists() or not sv_summary_file.exists():
        continue

    snps = 0
    indels = 0

    with open(stats_file) as handle:
        for line in handle:
            if not line.startswith("SN"):
                continue

            fields = line.rstrip().split("\t")

            if len(fields) < 4:
                continue

            label = fields[2].strip().rstrip(":")
            value = fields[3].strip()

            if label == "number of SNPs":
                snps = int(value)

            if label == "number of indels":
                indels = int(value)

    with open(sv_summary_file) as handle:
        sv_data = next(csv.DictReader(handle, delimiter="\t"))

    included = samples[:run_number]

    rows.append(
        {
            "run": f"run_{run_id}",
            "n_individuals": run_number,
            "n_haplotypes": run_number * 2,
            "individuals": ",".join(included),
            "reference_path": "Z0002_hap1_contig1",
            "snps": snps,
            "all_indels": indels,
            "sv_insertions_ge_50bp": int(
                sv_data["sv_insertions_ge_50bp"]
            ),
            "sv_deletions_ge_50bp": int(
                sv_data["sv_deletions_ge_50bp"]
            ),
            "inserted_bp_ge_50bp": int(
                sv_data["inserted_bp_ge_50bp"]
            ),
            "deleted_bp_ge_50bp": int(
                sv_data["deleted_bp_ge_50bp"]
            ),
            "largest_insertion_bp": int(
                sv_data["largest_insertion_bp"]
            ),
            "largest_deletion_bp": int(
                sv_data["largest_deletion_bp"]
            ),
        }
    )

df = pd.DataFrame(rows)

df.to_csv(
    results_dir / "pggb_region_sv_summary.tsv",
    sep="\t",
    index=False,
)

df.to_excel(
    results_dir / "pggb_region_sv_summary.xlsx",
    index=False,
)
