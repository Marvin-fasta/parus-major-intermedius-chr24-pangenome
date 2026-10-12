configfile: "config.yaml"


SAMPLE = config["sample"]
OUT = config["results_dir"]
THREADS = int(config["threads"])

INPUT_ASSEMBLY = config["input_assembly"]
ONT_READS = config["ont_reads"]
ILLUMINA_R1 = config["illumina_r1"]
ILLUMINA_R2 = config["illumina_r2"]

MEDAKA_MODEL = config["medaka_model"]

JAVA_HEAP = config["java_heap"]
PILON_CHUNKSIZE = config["pilon_chunksize"]
PILON_DOWNSAMPLE = config["pilon_downsample"]


rule all:
    input:
        f"{OUT}/medaka_pilon2/assembly.fasta",
        directory(f"{OUT}/hapdup_{SAMPLE}_polished")

# ------------------------------------------------------------
# Trim and filter Illumina reads
# ------------------------------------------------------------

rule fastp:
    input:
        r1=ILLUMINA_R1,
        r2=ILLUMINA_R2
    output:
        r1=f"{OUT}/illumina/{SAMPLE}_trimmed_R1.fq.gz",
        r2=f"{OUT}/illumina/{SAMPLE}_trimmed_R2.fq.gz",
        html=f"{OUT}/illumina/{SAMPLE}_fastp.html",
        json=f"{OUT}/illumina/{SAMPLE}_fastp.json"
    params:
        merged_r1=f"{OUT}/illumina/{SAMPLE}_merged_R1.fq.gz",
        merged_r2=f"{OUT}/illumina/{SAMPLE}_merged_R2.fq.gz"
    threads: THREADS
    shell:
        r"""
        mkdir -p {OUT}/illumina

        cat {input.r1} > {params.merged_r1}
        cat {input.r2} > {params.merged_r2}

        fastp \
            -i {params.merged_r1} \
            -I {params.merged_r2} \
            -o {output.r1} \
            -O {output.r2} \
            -w {threads} \
            --detect_adapter_for_pe \
            --html {output.html} \
            --json {output.json}
        """

# ------------------------------------------------------------
# Medaka
# Flye input assembly -> Medaka
# ------------------------------------------------------------

rule medaka:
    input:
        assembly=INPUT_ASSEMBLY,
        reads=ONT_READS
    output:
        assembly=f"{OUT}/medaka/consensus.fasta"
    params:
        outdir=f"{OUT}/medaka"
    threads: THREADS
    shell:
        r"""
        rm -rf {params.outdir}

        medaka_consensus \
            -i {input.reads} \
            -d {input.assembly} \
            -o {params.outdir} \
            -t {threads} \
            -m {MEDAKA_MODEL}
        """


# ------------------------------------------------------------
# Pilon round 1
# Flye -> Medaka -> Pilon
# ------------------------------------------------------------

rule pilon_round1:
    input:
        assembly=f"{OUT}/medaka/consensus.fasta",
        r1=f"{OUT}/illumina/{SAMPLE}_trimmed_R1.fq.gz",
        r2=f"{OUT}/illumina/{SAMPLE}_trimmed_R2.fq.gz"
    output:
        assembly=f"{OUT}/medaka_pilon1/assembly.fasta"
    params:
        outdir=f"{OUT}/medaka_pilon1",
        bam_full=f"{OUT}/medaka_pilon1/illumina_to_medaka.sorted.bam",
        bam_ds=f"{OUT}/medaka_pilon1/illumina_to_medaka.ds50.bam"
    threads: THREADS
    shell:
        r"""
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        bwa-mem2 index {input.assembly}

        bwa-mem2 mem -t {threads} \
            {input.assembly} \
            {input.r1} \
            {input.r2} | \
            samtools sort -@ {threads} \
            -o {params.bam_full}

        samtools view -@ {threads} \
            -b \
            -s {PILON_DOWNSAMPLE} \
            -o {params.bam_ds} \
            {params.bam_full}

        samtools index -@ {threads} {params.bam_ds}

        export JAVA_TOOL_OPTIONS="-Xmx{JAVA_HEAP}"

        pilon \
            --genome {input.assembly} \
            --frags {params.bam_ds} \
            --output assembly \
            --outdir {params.outdir} \
            --changes \
            --diploid \
            --fix snps,indels,gaps \
            --chunksize {PILON_CHUNKSIZE}
        """


# ------------------------------------------------------------
# Pilon round 2
# Flye -> Medaka -> Pilon -> Pilon
# ------------------------------------------------------------

rule pilon_round2:
    input:
        assembly=f"{OUT}/medaka_pilon1/assembly.fasta",
        r1=f"{OUT}/illumina/{SAMPLE}_trimmed_R1.fq.gz",
        r2=f"{OUT}/illumina/{SAMPLE}_trimmed_R2.fq.gz"
    output:
        assembly=f"{OUT}/medaka_pilon2/assembly.fasta"
    params:
        outdir=f"{OUT}/medaka_pilon2",
        bam_full=f"{OUT}/medaka_pilon2/illumina_to_pilon1.sorted.bam",
        bam_ds=f"{OUT}/medaka_pilon2/illumina_to_pilon1.ds50.bam"
    threads: THREADS
    shell:
        r"""
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        bwa-mem2 index {input.assembly}

        bwa-mem2 mem -t {threads} \
            {input.assembly} \
            {input.r1} \
            {input.r2} | \
            samtools sort -@ {threads} \
            -o {params.bam_full}

        samtools view -@ {threads} \
            -b \
            -s {PILON_DOWNSAMPLE} \
            -o {params.bam_ds} \
            {params.bam_full}

        samtools index -@ {threads} {params.bam_ds}

        export JAVA_TOOL_OPTIONS="-Xmx{JAVA_HEAP}"

        pilon \
            --genome {input.assembly} \
            --frags {params.bam_ds} \
            --output assembly \
            --outdir {params.outdir} \
            --changes \
            --diploid \
            --fix snps,indels,gaps \
            --chunksize {PILON_CHUNKSIZE}
        """
# ------------------------------------------------------------
# Map original ONT reads to the final polished assembly
# This BAM is required by HapDup
# ------------------------------------------------------------

rule hapdup_mapping:
    input:
        assembly=f"{OUT}/medaka_pilon2/assembly.fasta",
        reads=ONT_READS
    output:
        bam=f"{OUT}/hapdup_mapping_{SAMPLE}_polished.bam",
        bai=f"{OUT}/hapdup_mapping_{SAMPLE}_polished.bam.bai"
    threads: THREADS
    shell:
        r"""
        minimap2 -ax map-ont -t {threads} \
            {input.assembly} \
            {input.reads} | \
            samtools sort -@ {threads} \
                -o {output.bam}

        samtools index -@ {threads} {output.bam}
        """


# ------------------------------------------------------------
# HapDup
# Uses the polished assembly and the newly generated ONT BAM
# ------------------------------------------------------------

rule hapdup:
    input:
        assembly=f"{OUT}/medaka_pilon2/assembly.fasta",
        bam=f"{OUT}/hapdup_mapping_{SAMPLE}_polished.bam",
        bai=f"{OUT}/hapdup_mapping_{SAMPLE}_polished.bam.bai"
    output:
        outdir=directory(f"{OUT}/hapdup_{SAMPLE}_polished")
    params:
        results_dir=OUT
    threads: THREADS
    shell:
        r"""
        mkdir -p {params.results_dir}

        HD_DIR=$(realpath {params.results_dir})

        docker run --rm \
            -v "$HD_DIR:$HD_DIR" \
            -u "$(id -u):$(id -g)" \
            mkolmogo/hapdup:0.12 \
            hapdup \
            --assembly "$HD_DIR/medaka_pilon2/assembly.fasta" \
            --bam "$HD_DIR/hapdup_mapping_{SAMPLE}_polished.bam" \
            --out-dir "$HD_DIR/hapdup_{SAMPLE}_polished" \
            -t {threads} \
            --rtype ont
        """
