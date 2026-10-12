#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path


def reverse_complement(sequence):
    table = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )
    return sequence.translate(table)[::-1]


def read_fasta(path):
    sequences = {}
    name = None
    chunks = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)

                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        sequences[name] = "".join(chunks)

    return sequences


def parse_paf(path, min_mapq):
    hits = []

    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")

            if len(fields) < 12:
                continue

            query_name = fields[0]
            query_length = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            strand = fields[4]
            target_name = fields[5]
            target_length = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            matches = int(fields[9])
            alignment_length = int(fields[10])
            mapq = int(fields[11])

            tags = fields[12:]
            primary = any(tag == "tp:A:P" for tag in tags)

            if mapq < min_mapq:
                continue

            if not primary:
                continue

            query_covered = query_end - query_start

            hits.append(
                {
                    "query_name": query_name,
                    "query_length": query_length,
                    "query_start": query_start,
                    "query_end": query_end,
                    "query_covered": query_covered,
                    "strand": strand,
                    "target_name": target_name,
                    "target_length": target_length,
                    "target_start": target_start,
                    "target_end": target_end,
                    "matches": matches,
                    "alignment_length": alignment_length,
                    "mapq": mapq,
                }
            )

    if not hits:
        return None

    hits.sort(
        key=lambda x: (
            x["query_covered"],
            x["alignment_length"],
            x["matches"],
            x["mapq"],
        ),
        reverse=True,
    )

    return hits[0]


def write_fasta(path, header, sequence, width=80):
    with open(path, "w") as handle:
        handle.write(f">{header}\n")

        for i in range(0, len(sequence), width):
            handle.write(sequence[i:i + width] + "\n")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", required=True)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--hap", required=True)
    parser.add_argument("--flank", required=True, type=int)
    parser.add_argument("--min-mapq", required=True, type=int)
    parser.add_argument("--min-query-coverage", required=True, type=float)

    args = parser.parse_args()

    best = parse_paf(args.paf, args.min_mapq)

    if best is None:
        raise RuntimeError(
            "No primary PAF alignment passed the MAPQ threshold."
        )

    coverage = best["query_covered"] / best["query_length"]

    if coverage < args.min_query_coverage:
        raise RuntimeError(
            f"Best alignment covers only {coverage:.3f} of the query. "
            f"Required: {args.min_query_coverage:.3f}. "
            "The candidate region may be fragmented across phase blocks."
        )

    sequences = read_fasta(args.assembly)

    if best["target_name"] not in sequences:
        raise RuntimeError(
            f"Target contig {best['target_name']} not found in assembly FASTA."
        )

    contig = sequences[best["target_name"]]

    start_0 = max(0, best["target_start"] - args.flank)
    end_0 = min(len(contig), best["target_end"] + args.flank)

    extracted = contig[start_0:end_0]

    if best["strand"] == "-":
        extracted = reverse_complement(extracted)

    output_header = f"{args.sample}_hap{args.hap}_contig1"

    write_fasta(args.output, output_header, extracted)

    with open(args.metadata, "w") as handle:
        handle.write(
            "\t".join(
                [
                    "sample",
                    "haplotype",
                    "query_name",
                    "query_length",
                    "query_start_0based",
                    "query_end_0based",
                    "query_coverage_fraction",
                    "assembly_contig",
                    "contig_length",
                    "target_start_0based",
                    "target_end_0based",
                    "extraction_start_1based",
                    "extraction_end_1based",
                    "orientation",
                    "matches",
                    "alignment_length",
                    "mapq",
                ]
            )
            + "\n"
        )

        handle.write(
            "\t".join(
                [
                    args.sample,
                    args.hap,
                    best["query_name"],
                    str(best["query_length"]),
                    str(best["query_start"]),
                    str(best["query_end"]),
                    f"{coverage:.6f}",
                    best["target_name"],
                    str(best["target_length"]),
                    str(best["target_start"]),
                    str(best["target_end"]),
                    str(start_0 + 1),
                    str(end_0),
                    best["strand"],
                    str(best["matches"]),
                    str(best["alignment_length"]),
                    str(best["mapq"]),
                ]
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
