# Chromosome 24 Pangenome Workflow for *Parus major intermedius*

This repository contains Snakemake workflows developed for a bachelor thesis at TU Dortmund University. The workflows were used to construct and analyse a haplotype-resolved pangenome of chromosome 24 from 12 *Parus major intermedius* individuals.

## Overview

The workflow includes:

1. Long-read polishing with Medaka and optional Racon polishing
2. Short-read polishing with Pilon
3. Haplotype resolution with HapDup
4. Chromosome-wide pangenome graph construction with PGGB
5. Candidate-region pangenome graph construction
6. Variant deconstruction using VG and VCF-based analysis
7. Graph visualisation with ODGI

## Repository structure

```text
workflow/     Snakemake workflows
config/       Example configuration files
scripts/      Custom Python helper scripts
envs/         Software environment definitions
