#!/usr/bin/env bash
# This script annotates a VCF file using Ensembl Variant Effect Predictor (VEP).
# VEP is run in offline mode using the local cache and GRCh38 reference assembly.
# The output includes population allele frequency information, including gnomAD frequencies.

set -euo pipefail


# ___ INPUT AND OUTPUT ___

INPUT_VCF="/path/to/input.vcf.gz"
OUT_TSV="/path/to/output/VEP_gnomAD_annotation.tsv"


# ___ TOOL ___

VEP=vep


# ___ RUN VEP ___
# Annotate variants using the local VEP cache.
# The output includes gene symbols, canonical transcripts, HGVS notation,
# transcript biotype, variant class, population allele frequencies,
# gnomAD allele frequencies and maximum allele frequencies.

$VEP \
    -i "$INPUT_VCF" \
    -o "$OUT_TSV" \
    --cache \
    --offline \
    --assembly GRCh38 \
    --species homo_sapiens \
    --tab \
    --force_overwrite \
    --symbol \
    --canonical \
    --hgvs \
    --biotype \
    --variant_class \
    --af \
    --af_gnomad \
    --max_af \
    --pick_allele_gene \
    --allele_number

echo "Done: $OUT_TSV"
