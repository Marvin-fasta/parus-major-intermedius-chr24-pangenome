#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Read configuration.
# ------------------------------------------------------------

eval "$(
python - <<'PY'
import shlex
import yaml

with open("config.yaml", encoding="utf-8") as handle:
    c = yaml.safe_load(handle)

values = {
    "OUT": c["results_dir"],
    "THREADS": c["threads"],
    "REF_FASTA": c["reference_fasta"],
    "REF_CHROM": c["reference_chromosome"],
    "REGION_START": c["region_start"],
    "REGION_END": c["region_end"],
    "MIN_MAPQ": c["min_mapq"],
    "MIN_BLOCK_QUERY_LENGTH": c["min_block_query_length"],
    "PGGB_IDENTITY": c["pggb_identity"],
    "PGGB_SEGMENT_LENGTH": c["pggb_segment_length"],
}

for key, value in values.items():
    print(f"{key}={shlex.quote(str(value))}")
PY
)"

RUN_DIR="${OUT}/run_reference_1400000_1600000"
REFERENCE_DIR="${OUT}/reference"
MAPPING_DIR="${OUT}/mapping"
REGION_DIR="${OUT}/regions"

REFERENCE_REGION="${REFERENCE_DIR}/Parus_major_1400000_1600000.fasta"
INPUT_FASTA="${RUN_DIR}/input_haplotypes_plus_reference.fasta"

mkdir -p "${RUN_DIR}"
mkdir -p "${REFERENCE_DIR}"
mkdir -p "${MAPPING_DIR}"
mkdir -p "${REGION_DIR}"

echo "Reference interval:"
echo "${REF_CHROM}:${REGION_START}-${REGION_END}"

# ------------------------------------------------------------
# Extract exact local reference region.
# ------------------------------------------------------------

samtools faidx "${REF_FASTA}"

samtools faidx \
    "${REF_FASTA}" \
    "${REF_CHROM}:${REGION_START}-${REGION_END}" | \
awk '
    /^>/ {
        print ">Parus_major_reference_contig1"
        next
    }
    {
        print
    }
' > "${REFERENCE_REGION}"

test -s "${REFERENCE_REGION}"

echo "Reference region:"
seqkit stats "${REFERENCE_REGION}"

# ------------------------------------------------------------
# Create sample/haplotype/assembly list from config.yaml.
# ------------------------------------------------------------

python - <<'PY' > "${RUN_DIR}/haplotype_assemblies.tsv"
import yaml

with open("config.yaml", encoding="utf-8") as handle:
    c = yaml.safe_load(handle)

template = c["hapdup_template"]
special = c.get("hapdup_special", {})

print("sample\thaplotype\tassembly")

for sample in c["samples"]:
    for hap in (1, 2):
        if sample in special:
            assembly = special[sample][f"hap{hap}"]
        else:
            assembly = template.format(sample=sample, hap=hap)

        print(f"{sample}\t{hap}\t{assembly}")
PY

# ------------------------------------------------------------
# Map the exact 200 kb reference interval to every haplotype.
#
# The reference region is query; assembly is target.
# --secondary=no avoids redundant secondary alignments.
# ------------------------------------------------------------

tail -n +2 "${RUN_DIR}/haplotype_assemblies.tsv" | \
while IFS=$'\t' read -r sample hap assembly; do

    test -s "${assembly}"

    prefix="${sample}_hap${hap}"

    paf="${MAPPING_DIR}/${prefix}.paf"
    blocks="${MAPPING_DIR}/${prefix}.blocks.tsv"
    coverage="${MAPPING_DIR}/${prefix}.coverage.tsv"
    region_fasta="${REGION_DIR}/${prefix}.region.fasta"

    echo "Mapping ${prefix}"

    minimap2 \
        -x asm10 \
        --secondary=no \
        -t "${THREADS}" \
        "${assembly}" \
        "${REFERENCE_REGION}" \
        > "${paf}"

    python scripts/select_paf_blocks.py \
        --paf "${paf}" \
        --output "${blocks}" \
        --report "${coverage}" \
        --min-mapq "${MIN_MAPQ}" \
        --min-query-length "${MIN_BLOCK_QUERY_LENGTH}"

    n_blocks=$(awk 'END {print NR - 1}' "${blocks}")

    if [ "${n_blocks}" -lt 1 ]; then
        echo "ERROR: No usable mapped block for ${prefix}" >&2
        exit 1
    fi

    # Create a FASTA that may contain multiple records if the
    # reference interval is split across phaseblocks/contigs.
    : > "${region_fasta}"

    block_number=0

    tail -n +2 "${blocks}" | \
    while IFS=$'\t' read -r contig start end strand qstart qend mapq qlength; do

        block_number=$((block_number + 1))

        header="${prefix}_block${block_number}"

        if [ "${strand}" = "-" ]; then
            samtools faidx \
                "${assembly}" \
                "${contig}:${start}-${end}" | \
            seqkit seq -r -p | \
            awk -v header=">${header}" '
                NR == 1 {
                    print header
                    next
                }
                {
                    print
                }
            ' >> "${region_fasta}"
        else
            samtools faidx \
                "${assembly}" \
                "${contig}:${start}-${end}" | \
            awk -v header=">${header}" '
                NR == 1 {
                    print header
                    next
                }
                {
                    print
                }
            ' >> "${region_fasta}"
        fi
    done

    test -s "${region_fasta}"

done

# ------------------------------------------------------------
# Store block-coverage summary for quality control.
# ------------------------------------------------------------

{
    echo -e "sample\thaplotype\tquery_length_bp\tcovered_query_bp\tcoverage_percent\tblocks"

    tail -n +2 "${RUN_DIR}/haplotype_assemblies.tsv" | \
    while IFS=$'\t' read -r sample hap assembly; do
        coverage="${MAPPING_DIR}/${sample}_hap${hap}.coverage.tsv"

        tail -n +2 "${coverage}" | \
        awk -v sample="${sample}" -v hap="${hap}" \
            'BEGIN {OFS="\t"} {print sample, hap, $0}'
    done

} > "${RUN_DIR}/mapping_coverage_summary.tsv"

echo
echo "Mapping coverage summary:"
column -t -s $'\t' "${RUN_DIR}/mapping_coverage_summary.tsv"

# ------------------------------------------------------------
# Build PGGB input:
# reference first, then all haplotype blocks.
# ------------------------------------------------------------

cat "${REFERENCE_REGION}" > "${INPUT_FASTA}"

tail -n +2 "${RUN_DIR}/haplotype_assemblies.tsv" | \
while IFS=$'\t' read -r sample hap assembly; do
    cat "${REGION_DIR}/${sample}_hap${hap}.region.fasta"
done >> "${INPUT_FASTA}"

samtools faidx "${INPUT_FASTA}"

grep '^>' "${INPUT_FASTA}" \
    > "${RUN_DIR}/input_sequence_headers.txt"

N_FASTA_RECORDS=$(grep -c '^>' "${INPUT_FASTA}")

echo "FASTA records supplied to PGGB: ${N_FASTA_RECORDS}"

grep -Fx ">Parus_major_reference_contig1" \
    "${INPUT_FASTA}"

# 12 diploid individuals = 24 haplotypes, plus reference = 25.
N_HAPLOTYPES=25

# ------------------------------------------------------------
# Run PGGB.
# ------------------------------------------------------------

RAW_DIR="${RUN_DIR}/pggb_raw"
LOG="${RUN_DIR}/pggb.log"

rm -rf "${RAW_DIR}"

pggb \
    -i "${INPUT_FASTA}" \
    -o "${RAW_DIR}" \
    -n "${N_HAPLOTYPES}" \
    -p "${PGGB_IDENTITY}" \
    -s "${PGGB_SEGMENT_LENGTH}" \
    -t "${THREADS}" \
    > "${LOG}" 2>&1

GFA=$(find "${RAW_DIR}" -type f -name "*.smooth.final.gfa" | head -n 1)

test -n "${GFA}"
test -s "${GFA}"

cp "${GFA}" "${RUN_DIR}/graph.smooth.final.gfa"

odgi build \
    -g "${RUN_DIR}/graph.smooth.final.gfa" \
    -o "${RUN_DIR}/graph.smooth.final.og"

# ------------------------------------------------------------
# VG conversion, XG index and VCF relative to Parus major reference.
# ------------------------------------------------------------

vg convert \
    -g "${RUN_DIR}/graph.smooth.final.gfa" \
    > "${RUN_DIR}/graph.vg"

vg index \
    -x "${RUN_DIR}/graph.xg" \
    "${RUN_DIR}/graph.vg"

vg deconstruct \
    -p "Parus_major_reference_contig1" \
    "${RUN_DIR}/graph.vg" \
    > "${RUN_DIR}/variants.vcf" \
    2> "${RUN_DIR}/vg_deconstruct.log"

bcftools sort \
    "${RUN_DIR}/variants.vcf" \
    -Oz \
    -o "${RUN_DIR}/variants.sorted.vcf.gz"

tabix -f -p vcf "${RUN_DIR}/variants.sorted.vcf.gz"

bcftools stats \
    "${RUN_DIR}/variants.sorted.vcf.gz" \
    > "${RUN_DIR}/variants.stats.txt"

# ------------------------------------------------------------
# Variant tables for later SNP/Indel plotting.
# ------------------------------------------------------------

python scripts/extract_all_variants.py \
    --vcf "${RUN_DIR}/variants.sorted.vcf.gz" \
    --output "${RUN_DIR}/all_variants.tsv"

python scripts/make_snp_windows.py \
    --variants "${RUN_DIR}/all_variants.tsv" \
    --output "${RUN_DIR}/snp_density_1kb.tsv" \
    --reference-contig "Parus_major_reference_contig1" \
    --region-length 200001 \
    --window-size 1000

python scripts/make_snp_windows.py \
    --variants "${RUN_DIR}/all_variants.tsv" \
    --output "${RUN_DIR}/snp_density_5kb.tsv" \
    --reference-contig "Parus_major_reference_contig1" \
    --region-length 200001 \
    --window-size 5000

# ------------------------------------------------------------
# Graph statistics and visualisation.
# ------------------------------------------------------------

odgi stats \
    -i "${RUN_DIR}/graph.smooth.final.og" \
    > "${RUN_DIR}/odgi_stats.txt"

odgi viz \
    -i "${RUN_DIR}/graph.smooth.final.og" \
    -o "${RUN_DIR}/odgi_viz.png"

touch "${RUN_DIR}/run_complete.done"

echo
echo "============================================"
echo "PGGB analysis successfully completed."
echo "Results: ${RUN_DIR}"
echo "============================================"
