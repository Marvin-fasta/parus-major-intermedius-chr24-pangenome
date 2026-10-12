configfile: "config.yaml"

from snakemake.io import directory


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

MERQURY_K = config["merqury_k"]


RACON_INPUTS = {
    "racon1": INPUT_ASSEMBLY,
    "racon2": f"{OUT}/racon1/assembly.fasta"
}


MEDAKA_INPUTS = {
    "medaka": INPUT_ASSEMBLY,
    "racon2_medaka": f"{OUT}/racon2/assembly.fasta"
}


PILON_INPUTS = {
    "pilon1": INPUT_ASSEMBLY,
    "pilon2": f"{OUT}/pilon1/assembly.fasta",
    "medaka_pilon1": f"{OUT}/medaka/consensus.fasta",
    "medaka_pilon2": f"{OUT}/medaka_pilon1/assembly.fasta",
    "racon2_medaka_pilon1": f"{OUT}/racon2_medaka/consensus.fasta",
    "racon2_medaka_pilon2": f"{OUT}/racon2_medaka_pilon1/assembly.fasta"
}


MERQURY_ASSEMBLIES = {
    "flye": INPUT_ASSEMBLY,
    "flye_racon1": f"{OUT}/racon1/assembly.fasta",
    "flye_racon2": f"{OUT}/racon2/assembly.fasta",
    "flye_racon2_medaka": f"{OUT}/racon2_medaka/consensus.fasta",
    "flye_medaka": f"{OUT}/medaka/consensus.fasta",
    "flye_medaka_pilon1": f"{OUT}/medaka_pilon1/assembly.fasta",
    "flye_medaka_pilon2": f"{OUT}/medaka_pilon2/assembly.fasta",
    "flye_racon2_medaka_pilon1": f"{OUT}/racon2_medaka_pilon1/assembly.fasta",
    "flye_racon2_medaka_pilon2": f"{OUT}/racon2_medaka_pilon2/assembly.fasta",
    "flye_pilon1": f"{OUT}/pilon1/assembly.fasta",
    "flye_pilon2": f"{OUT}/pilon2/assembly.fasta"
}


rule all:
    input:
        expand(
            f"{OUT}/merqury/{{variant}}/.done",
            variant=MERQURY_ASSEMBLIES.keys()
        )


rule fastp:
    input:
        r1=ILLUMINA_R1,
        r2=ILLUMINA_R2
    output:
        r1=f"{OUT}/illumina/{SAMPLE}_trimmed_R1.fq.gz",
        r2=f"{OUT}/illumina/{SAMPLE}_trimmed_R2.fq.gz",
        html=f"{OUT}/illumina/{SAMPLE}_fastp.html",
        json=f"{OUT}/illumina/{SAMPLE}_fastp.json"
    threads: THREADS
    shell:
        r"""
        mkdir -p {OUT}/illumina

        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            -w {threads} \
            --detect_adapter_for_pe \
            --html {output.html} \
            --json {output.json}
        """


rule racon:
    input:
        assembly=lambda wc: RACON_INPUTS[wc.route],
        reads=ONT_READS
    output:
        assembly=f"{OUT}/{{route}}/assembly.fasta"
    params:
        outdir=lambda wc: f"{OUT}/{wc.route}",
        paf=lambda wc: f"{OUT}/{wc.route}/ont_alignment.paf",
        paf_filtered=lambda wc: f"{OUT}/{wc.route}/ont_alignment.primary.q20.paf"
    threads: THREADS
    wildcard_constraints:
        route="racon1|racon2"
    shell:
        r"""
        mkdir -p {params.outdir}

        minimap2 -x map-ont -t {threads} \
            {input.assembly} \
            {input.reads} \
            > {params.paf}

        awk '$12 >= 20 && $0 ~ /\ttp:A:P(\t|$)/' \
            {params.paf} \
            > {params.paf_filtered}

        racon -t {threads} \
            {input.reads} \
            {params.paf_filtered} \
            {input.assembly} \
            > {output.assembly}
        """


rule medaka:
    input:
        assembly=lambda wc: MEDAKA_INPUTS[wc.route],
        reads=ONT_READS
    output:
        assembly=f"{OUT}/{{route}}/consensus.fasta"
    params:
        outdir=lambda wc: f"{OUT}/{wc.route}"
    threads: THREADS
    wildcard_constraints:
        route="medaka|racon2_medaka"
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


rule pilon:
    input:
        assembly=lambda wc: PILON_INPUTS[wc.route],
        r1=f"{OUT}/illumina/{SAMPLE}_trimmed_R1.fq.gz",
        r2=f"{OUT}/illumina/{SAMPLE}_trimmed_R2.fq.gz"
    output:
        assembly=f"{OUT}/{{route}}/assembly.fasta"
    params:
        outdir=lambda wc: f"{OUT}/{wc.route}",
        bam_full=lambda wc: f"{OUT}/{wc.route}/illumina_to_assembly.sorted.bam",
        bam_ds=lambda wc: f"{OUT}/{wc.route}/illumina_to_assembly.ds50.bam"
    threads: THREADS
    wildcard_constraints:
        route="pilon1|pilon2|medaka_pilon1|medaka_pilon2|racon2_medaka_pilon1|racon2_medaka_pilon2"
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


rule meryl_count:
    input:
        r1=f"{OUT}/illumina/{SAMPLE}_trimmed_R1.fq.gz",
        r2=f"{OUT}/illumina/{SAMPLE}_trimmed_R2.fq.gz"
    output:
        db=directory(f"{OUT}/merqury/{SAMPLE}_illumina_k{MERQURY_K}.meryl")
    params:
        k=MERQURY_K
    shell:
        r"""
        mkdir -p {OUT}/merqury
        rm -rf {output.db}

        meryl k={params.k} count \
            {input.r1} \
            {input.r2} \
            output {output.db}
        """


rule merqury:
    input:
        db=f"{OUT}/merqury/{SAMPLE}_illumina_k{MERQURY_K}.meryl",
        assembly=lambda wc: MERQURY_ASSEMBLIES[wc.variant]
    output:
        done=f"{OUT}/merqury/{{variant}}/.done"
    params:
        outdir=lambda wc: f"{OUT}/merqury/{wc.variant}"
    wildcard_constraints:
        variant="flye|flye_racon1|flye_racon2|flye_racon2_medaka|flye_medaka|flye_medaka_pilon1|flye_medaka_pilon2|flye_racon2_medaka_pilon1|flye_racon2_medaka_pilon2|flye_pilon1|flye_pilon2"
    shell:
        r"""
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        DB=$(realpath {input.db})
        ASSEMBLY=$(realpath {input.assembly})

        (
            cd {params.outdir}

            merqury.sh \
                "$DB" \
                "$ASSEMBLY" \
                "{wildcards.variant}"
        )

        touch {output.done}
        """
