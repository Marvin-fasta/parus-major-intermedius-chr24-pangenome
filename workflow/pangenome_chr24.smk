configfile: "config.yaml"

from pathlib import Path


OUT = config["results_dir"]
THREADS = int(config["threads"])

SAMPLES = config["samples"]

PGGB_IDENTITY = int(config["pggb_identity"])
PGGB_SEGMENT_LENGTH = int(config["pggb_segment_length"])

SV_MIN_SIZE = int(config["sv_min_size"])

REFERENCE_SAMPLE = config["reference_sample"]
REFERENCE_HAPLOTYPE = int(config["reference_haplotype"])

HAPDUP_TEMPLATE = config["hapdup_template"]
HAPDUP_SPECIAL = config.get("hapdup_special", {})

RUN_IDS = [f"{i:02d}" for i in range(1, len(SAMPLES) + 1)]


def hapdup_fasta(sample, haplotype):
    if sample in HAPDUP_SPECIAL:
        return HAPDUP_SPECIAL[sample][f"hap{haplotype}"]

    return HAPDUP_TEMPLATE.format(
        sample=sample,
        hap=haplotype
    )


def samples_for_run(run_id):
    return SAMPLES[:int(run_id)]


def renamed_inputs_for_run(wildcards):
    sample_list = samples_for_run(wildcards.run)

    files = []

    for sample in sample_list:
        files.append(f"{OUT}/renamed/{sample}_hap1.fasta")
        files.append(f"{OUT}/renamed/{sample}_hap2.fasta")

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
        f"{OUT}/pggb_sv_summary.xlsx",
        f"{OUT}/pggb_sv_summary.tsv"


# ------------------------------------------------------------
# Rename all contigs from each HapDup phased FASTA.
#
# Example:
# >Z0072_contig_13_hap1
# >Z0072_hap1_contig1
#
# Every input FASTA receives unique headers.
# ------------------------------------------------------------

rule rename_haplotype_fasta:
    input:
        fasta=lambda wc: hapdup_fasta(
            wc.sample,
            int(wc.hap)
        )
    output:
        fasta=f"{OUT}/renamed/{{sample}}_hap{{hap}}.fasta"
    shell:
        r"""
        mkdir -p {OUT}/renamed

        seqkit replace \
            -p "^.+$" \
            -r "{wildcards.sample}_hap{wildcards.hap}_contig{{nr}}" \
            {input.fasta} \
            > {output.fasta}
        """


# ------------------------------------------------------------
# Create cumulative FASTA input for each PGGB run.
#
# run_01 = Z0002
# run_02 = Z0002 + Z0016
# ...
# run_12 = all samples
# ------------------------------------------------------------

rule merge_cumulative_haplotypes:
    input:
        renamed_inputs_for_run
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
# Run PGGB.
#
# Each PGGB run depends on completion of the previous full run.
# Therefore run_01 -> run_02 -> ... -> run_12 is sequential.
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
            -t {threads}

        GFA=$(find {params.raw_dir} -type f -name "*.smooth.final.gfa" | head -n 1)
        OG=$(find {params.raw_dir} -type f -name "*.smooth.final.og" | head -n 1)

        test -n "$GFA"
        test -n "$OG"

        cp "$GFA" {output.gfa}
        cp "$OG" {output.og}

        touch {output.done}
        """


# ------------------------------------------------------------
# Convert PGGB GFA to VG and create XG index.
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
# Get all reference phase blocks from Z0002 hap1.
#
# This covers all Z0002_hap1_contig1, contig2, contig3 ...
# instead of only one contig.
# ------------------------------------------------------------

rule reference_paths:
    input:
        xg=f"{OUT}/runs/run_{{run}}/graph.xg"
    output:
        paths=f"{OUT}/runs/run_{{run}}/reference_paths.txt"
    shell:
        r"""
        vg paths -L -x {input.xg} | \
            grep "^Z0002_hap1_contig" \
            > {output.paths}

        test -s {output.paths}
        """


# ------------------------------------------------------------
# Deconstruct graph relative to all Z0002 hap1 phase blocks.
#
# A VCF is generated separately for each reference block.
# Then all VCFs are concatenated into one VCF.
# ------------------------------------------------------------

rule vg_deconstruct:
    input:
        vg=f"{OUT}/runs/run_{{run}}/graph.vg",
        paths=f"{OUT}/runs/run_{{run}}/reference_paths.txt"
    output:
        vcf=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz",
        tbi=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz.tbi"
    params:
        block_dir=lambda wc: f"{OUT}/runs/run_{wc.run}/deconstruct_blocks",
        sites_dir=lambda wc: f"{OUT}/runs/run_{wc.run}/deconstruct_sites_only"
    log:
        f"{OUT}/runs/run_{{run}}/vg_deconstruct.log"
    shell:
        r"""
        rm -rf {params.block_dir}
        rm -rf {params.sites_dir}

        mkdir -p {params.block_dir}
        mkdir -p {params.sites_dir}

        while read -r REF_PATH
        do
            SAFE_NAME=$(echo "$REF_PATH" | tr '/: ' '___')

            # Original VCF with genotype/sample information.
            vg deconstruct \
                -p "$REF_PATH" \
                {input.vg} \
                > {params.block_dir}/"$SAFE_NAME".vcf

            bcftools sort \
                {params.block_dir}/"$SAFE_NAME".vcf \
                -Oz \
                -o {params.block_dir}/"$SAFE_NAME".vcf.gz

            tabix -f -p vcf \
                {params.block_dir}/"$SAFE_NAME".vcf.gz

            # Remove genotype/sample columns for cross-contig summary.
            bcftools view \
                -G \
                -Oz \
                -o {params.sites_dir}/"$SAFE_NAME".sites.vcf.gz \
                {params.block_dir}/"$SAFE_NAME".vcf.gz

            tabix -f -p vcf \
                {params.sites_dir}/"$SAFE_NAME".sites.vcf.gz

        done < {input.paths}

        VCF_FILES=$(find {params.sites_dir} \
            -type f \
            -name "*.sites.vcf.gz" | sort)

        test -n "$VCF_FILES"

        bcftools concat \
            -a \
            --allow-overlaps \
            -Oz \
            -o {output.vcf} \
            $VCF_FILES \
            > {log} 2>&1

        tabix -f -p vcf {output.vcf}
        """

# ------------------------------------------------------------
# bcftools stats and plots.
# ------------------------------------------------------------

rule bcftools_stats:
    input:
        vcf=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz",
        tbi=f"{OUT}/runs/run_{{run}}/variants.sorted.vcf.gz.tbi"
    output:
        stats=f"{OUT}/runs/run_{{run}}/variants.stats.txt",
        stats_done=f"{OUT}/runs/run_{{run}}/bcftools_stats.done"
    params:
        plot_dir=lambda wc: f"{OUT}/runs/run_{wc.run}/vcfstats_plots"
    shell:
        r"""
        bcftools stats {input.vcf} > {output.stats}

        rm -rf {params.plot_dir}

        plot-vcfstats \
            -p {params.plot_dir} \
            {output.stats}

        touch {output.stats_done}
        """


# ------------------------------------------------------------
# Extract SV candidates:
# insertions/deletions with absolute size >= SV_MIN_SIZE.
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
# ODGI graph statistics and graph visualization.
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
# Marker: the complete run is finished.
#
# The next PGGB run waits for this marker.
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
# Final Excel and TSV summary over all cumulative runs.
# ------------------------------------------------------------

rule excel_summary:
    input:
        complete=expand(
            f"{OUT}/runs/run_{{run}}/run_complete.done",
            run=RUN_IDS
        )
    output:
        xlsx=f"{OUT}/pggb_sv_summary.xlsx",
        tsv=f"{OUT}/pggb_sv_summary.tsv"
    params:
        samples=",".join(SAMPLES)
    shell:
        r"""
        python scripts/make_excel_summary.py \
            {OUT} \
            "{params.samples}"
        """
