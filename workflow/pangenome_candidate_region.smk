configfile: "config.yaml"


OUT = config["results_dir"]
THREADS = int(config["threads"])

SAMPLES = config["samples"]

REFERENCE_FASTA = config["reference_fasta"]
REFERENCE_CHROM = config["reference_chromosome"]
REGION_START = int(config["region_start"])
REGION_END = int(config["region_end"])

FLANK_BP = int(config["flank_bp"])
MIN_QUERY_COVERAGE = float(config["min_query_coverage"])
MIN_MAPQ = int(config["min_mapq"])

PGGB_IDENTITY = int(config["pggb_identity"])
PGGB_SEGMENT_LENGTH = int(config["pggb_segment_length"])

SV_MIN_SIZE = int(config["sv_min_size"])

HAPDUP_TEMPLATE = config["hapdup_template"]
HAPDUP_SPECIAL = config.get("hapdup_special", {})

RUN_IDS = [f"{i:02d}" for i in range(1, len(SAMPLES) + 1)]
MANUAL_REGION_DIR = config["manual_region_dir"]
)

MANUAL_REGIONS = {
    ("Z0064", "1"): f"{MANUAL_REGION_DIR}/Z0064_hap1_candidate_region.fasta",
    ("Z0064", "2"): f"{MANUAL_REGION_DIR}/Z0064_hap2_candidate_region.fasta",
    ("Z0043", "1"): f"{MANUAL_REGION_DIR}/Z0043_hap1_candidate_region.fasta",
    ("Z0043", "2"): f"{MANUAL_REGION_DIR}/Z0043_hap2_candidate_region.fasta",
    ("Z0049", "1"): f"{MANUAL_REGION_DIR}/Z0049_hap1_candidate_region.fasta",
}

def hapdup_fasta(sample, haplotype):
    if sample in HAPDUP_SPECIAL:
        return HAPDUP_SPECIAL[sample][f"hap{haplotype}"]

    return HAPDUP_TEMPLATE.format(
        sample=sample,
        hap=haplotype,
    )


def samples_for_run(run_id):
    return SAMPLES[:int(run_id)]


def regional_fastas_for_run(wildcards):
    files = []

    for sample in samples_for_run(wildcards.run):
        files.append(f"{OUT}/regions/{sample}_hap1.region.fasta")
        files.append(f"{OUT}/regions/{sample}_hap2.region.fasta")

    return files


def previous_run_complete(wildcards):
    current = int(wildcards.run)

    if current == 1:
        return []

    previous = f"{current - 1:02d}"

    return [
        f"{OUT}/runs/run_{previous}/run_complete.done"
    ]


rule all:
    input:
        f"{OUT}/pggb_region_sv_summary.xlsx",
        f"{OUT}/pggb_region_sv_summary.tsv"


# ------------------------------------------------------------
# Index Parus major reference genome.
# ------------------------------------------------------------

rule index_reference:
    input:
        fasta=REFERENCE_FASTA
    output:
        fai=REFERENCE_FASTA + ".fai"
    shell:
        r"""
        samtools faidx {input.fasta}
        """


# ------------------------------------------------------------
# Extract exact 34,925 bp candidate region from Parus major.
# ------------------------------------------------------------

rule extract_reference_region:
    input:
        fasta=REFERENCE_FASTA,
        fai=REFERENCE_FASTA + ".fai"
    output:
        fasta=f"{OUT}/reference/Parus_major_candidate_exact.fasta"
    params:
        region=f"{REFERENCE_CHROM}:{REGION_START}-{REGION_END}"
    shell:
        r"""
        mkdir -p {OUT}/reference

        samtools faidx \
            {input.fasta} \
            {params.region} \
            > {output.fasta}
        """


# ------------------------------------------------------------
# Map exact reference candidate interval to every haplotype.
# Assembly is the target, Parus major interval is the query.
# ------------------------------------------------------------

rule map_reference_region:
    input:
        assembly=lambda wc: hapdup_fasta(
            wc.sample,
            int(wc.hap)
        ),
        region=f"{OUT}/reference/Parus_major_candidate_exact.fasta"
    output:
        paf=f"{OUT}/mapping/{{sample}}_hap{{hap}}.paf"
    threads: THREADS
    shell:
        r"""
        mkdir -p {OUT}/mapping

        minimap2 \
            -x asm10 \
            -t {threads} \
            {input.assembly} \
            {input.region} \
            > {output.paf}
        """


# ------------------------------------------------------------
# Extract best mapped region plus/minus 10 kb.
# Reverse-complement output automatically if PAF orientation is '-'.
# ------------------------------------------------------------

rule extract_haplotype_region:
    input:
        assembly=lambda wc: hapdup_fasta(
            wc.sample,
            int(wc.hap)
        ),
        paf=f"{OUT}/mapping/{{sample}}_hap{{hap}}.paf",
        manual=lambda wc: (
            [MANUAL_REGIONS[(wc.sample, wc.hap)]]
            if (wc.sample, wc.hap) in MANUAL_REGIONS
            else []
        )
    output:
        fasta=f"{OUT}/regions/{{sample}}_hap{{hap}}.region.fasta",
        metadata=f"{OUT}/regions/{{sample}}_hap{{hap}}.region.tsv"
    params:
        script="scripts/extract_region_from_paf.py"
    shell:
        r"""
        mkdir -p {OUT}/regions

        if [ -n "{input.manual}" ]; then

            # The manually extracted FASTA can contain multiple phase blocks.
            # Rename all records consistently for PGGB.
            seqkit replace \
                -p "^.+$" \
                -r "{wildcards.sample}_hap{wildcards.hap}_contig{{nr}}" \
                {input.manual} \
                > {output.fasta}

            printf "sample\thaplotype\textraction_method\tmanual_file\n" \
                > {output.metadata}

            printf "{wildcards.sample}\t{wildcards.hap}\tmanual_phaseblocks\t{input.manual}\n" \
                >> {output.metadata}

        else

            python {params.script} \
                --assembly {input.assembly} \
                --paf {input.paf} \
                --output {output.fasta} \
                --metadata {output.metadata} \
                --sample {wildcards.sample} \
                --hap {wildcards.hap} \
                --flank {FLANK_BP} \
                --min-mapq {MIN_MAPQ} \
                --min-query-coverage {MIN_QUERY_COVERAGE}

        fi
        """
# ------------------------------------------------------------
# Merge the regional haplotypes for every cumulative run.
# ------------------------------------------------------------

rule merge_cumulative_regions:
    input:
        regional_fastas_for_run
    output:
        fasta=f"{OUT}/runs/run_{{run}}/input_haplotypes.fasta"
    params:
        sample_names=lambda wc: ",".join(samples_for_run(wc.run))
    shell:
        r"""
        mkdir -p {OUT}/runs/run_{wildcards.run}

        cat {input} > {output.fasta}

        echo "{params.sample_names}" \
            > {OUT}/runs/run_{wildcards.run}/included_samples.txt
        """


rule index_pggb_input:
    input:
        fasta=f"{OUT}/runs/run_{{run}}/input_haplotypes.fasta"
    output:
        fai=f"{OUT}/runs/run_{{run}}/input_haplotypes.fasta.fai"
    shell:
        r"""
        samtools faidx {input.fasta}
        """


# ------------------------------------------------------------
# PGGB is sequential: run 02 starts only after run 01 is complete.
# ------------------------------------------------------------

rule pggb:
    input:
        fasta=f"{OUT}/runs/run_{{run}}/input_haplotypes.fasta",
        fai=f"{OUT}/runs/run_{{run}}/input_haplotypes.fasta.fai",
        previous=previous_run_complete
    output:
        gfa=f"{OUT}/runs/run_{{run}}/graph.smooth.final.gfa",
        og=f"{OUT}/runs/run_{{run}}/graph.smooth.final.og",
        done=f"{OUT}/runs/run_{{run}}/pggb.done"
    params:
        raw_dir=lambda wc: f"{OUT}/runs/run_{wc.run}/pggb_raw",
        n_haplotypes=lambda wc: 2 * int(wc.run)
    threads: THREADS
    log:
        f"{OUT}/runs/run_{{run}}/pggb.log"
    shell:
        r"""
        rm -rf {params.raw_dir}
        mkdir -p {OUT}/runs/run_{wildcards.run}

        pggb \
            -i {input.fasta} \
            -o {params.raw_dir} \
            -n {params.n_haplotypes} \
            -p {PGGB_IDENTITY} \
            -s {PGGB_SEGMENT_LENGTH} \
            -t {threads} \
            > {log} 2>&1

        GFA=$(find {params.raw_dir} -type f -name "*.smooth.final.gfa" | head -n 1)

        test -n "$GFA"
        test -s "$GFA"

        cp "$GFA" {output.gfa}

        odgi build \
            -g {output.gfa} \
            -o {output.og}

        touch {output.done}
        """


# ------------------------------------------------------------
# Convert final graph to VG and build XG index.
# ------------------------------------------------------------

rule vg_convert_and_index:
    input:
        gfa=f"{OUT}/runs/run_{{run}}/graph.smooth.final.gfa"
    output:
        vg=f"{OUT}/runs/run_{{run}}/graph.vg",
        xg=f"{OUT}/runs/run_{{run}}/graph.xg"
    shell:
        r"""
        vg convert -g {input.gfa} > {output.vg}

        vg index \
            -x {output.xg} \
            {output.vg}
        """


# ------------------------------------------------------------
# All regional sequences have one contig each.
# Z0002 hap1 is always the VCF reference path.
# ------------------------------------------------------------

rule vg_deconstruct:
    input:
        vg=f"{OUT}/runs/run_{{run}}/graph.vg"
    output:
        vcf=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz",
        tbi=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz.tbi"
    params:
        reference_path="Z0002_hap1_contig1"
    log:
        f"{OUT}/runs/run_{{run}}/vg_deconstruct.log"
    shell:
        r"""
        vg deconstruct \
            -p {params.reference_path} \
            {input.vg} \
            > {OUT}/runs/run_{wildcards.run}/variants.vcf \
            2> {log}

        bcftools sort \
            {OUT}/runs/run_{wildcards.run}/variants.vcf \
            -Oz \
            -o {output.vcf}

        tabix -f -p vcf {output.vcf}
        """


# ------------------------------------------------------------
# VCF statistics. No plot-vcfstats here, to avoid PDF failures.
# ------------------------------------------------------------

rule bcftools_stats:
    input:
        vcf=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz",
        tbi=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz.tbi"
    output:
        stats=f"{OUT}/runs/run_{{run}}/variants.stats.txt"
    shell:
        r"""
        bcftools stats {input.vcf} > {output.stats}
        """


# ------------------------------------------------------------
# Extract insertion and deletion SVs >= 50 bp.
# ------------------------------------------------------------

rule sv_table:
    input:
        vcf=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz",
        tbi=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz.tbi"
    output:
        table=f"{OUT}/runs/run_{{run}}/sv_ge_50bp.tsv",
        summary=f"{OUT}/runs/run_{{run}}/sv_summary.tsv"
    shell:
        r"""
        echo -e "CHROM\tPOS\tREF_length\tALT_length\tSV_length\tSV_type" \
            > {output.table}

        bcftools norm -m -any -Ou {input.vcf} | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' | \
        awk -F'\t' -v min_size={SV_MIN_SIZE} '
        BEGIN {{
            OFS="\t"
        }}
        {{
            ref = $3
            alt = $4

            if (alt == "." || alt == "*" || alt ~ /^</) {{
                next
            }}

            ref_len = length(ref)
            alt_len = length(alt)
            diff = alt_len - ref_len

            if (diff >= min_size) {{
                print $1, $2, ref_len, alt_len, diff, "INS"
            }}
            else if (diff <= -min_size) {{
                print $1, $2, ref_len, alt_len, diff, "DEL"
            }}
        }}' \
        >> {output.table}

        awk -F'\t' '
        BEGIN {{
            OFS="\t"
            ins = 0
            del = 0
            ins_bp = 0
            del_bp = 0
            max_ins = 0
            max_del = 0
        }}
        NR > 1 {{
            if ($6 == "INS") {{
                ins++
                ins_bp += $5

                if ($5 > max_ins) {{
                    max_ins = $5
                }}
            }}
            else if ($6 == "DEL") {{
                del++
                abs_len = -$5
                del_bp += abs_len

                if (abs_len > max_del) {{
                    max_del = abs_len
                }}
            }}
        }}
        END {{
            print "sv_insertions_ge_50bp", \
                  "sv_deletions_ge_50bp", \
                  "inserted_bp_ge_50bp", \
                  "deleted_bp_ge_50bp", \
                  "largest_insertion_bp", \
                  "largest_deletion_bp"

            print ins, del, ins_bp, del_bp, max_ins, max_del
        }}' \
        {output.table} \
        > {output.summary}
        """


# ------------------------------------------------------------
# ODGI graph statistics and visualization.
# ------------------------------------------------------------

rule odgi_outputs:
    input:
        og=f"{OUT}/runs/run_{{run}}/graph.smooth.final.og"
    output:
        stats=f"{OUT}/runs/run_{{run}}/odgi_stats.txt",
        viz=f"{OUT}/runs/run_{{run}}/odgi_viz.png"
    shell:
        r"""
        odgi stats \
            -i {input.og} \
            > {output.stats}

        odgi viz \
            -i {input.og} \
            -o {output.viz}
        """


# ------------------------------------------------------------
# Completion marker. The next cumulative PGGB run waits here.
# ------------------------------------------------------------

rule run_complete:
    input:
        stats=f"{OUT}/runs/run_{{run}}/variants.stats.txt",
        sv_summary=f"{OUT}/runs/run_{{run}}/sv_summary.tsv",
        odgi_stats=f"{OUT}/runs/run_{{run}}/odgi_stats.txt",
        odgi_viz=f"{OUT}/runs/run_{{run}}/odgi_viz.png"
    output:
        done=f"{OUT}/runs/run_{{run}}/run_complete.done"
    shell:
        r"""
        touch {output.done}
        """


# ------------------------------------------------------------
# Final Excel and TSV summary.
# ------------------------------------------------------------

rule excel_summary:
    input:
        complete=expand(
            f"{OUT}/runs/run_{{run}}/run_complete.done",
            run=RUN_IDS
        )
    output:
        xlsx=f"{OUT}/pggb_region_sv_summary.xlsx",
        tsv=f"{OUT}/pggb_region_sv_summary.tsv"
    params:
        samples=",".join(SAMPLES)
    shell:
        r"""
        python scripts/make_excel_summary.py \
            {OUT} \
            "{params.samples}"
        """
