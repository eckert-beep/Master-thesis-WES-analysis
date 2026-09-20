#!/usr/bin/env bash
# This script annotates a VCF file using Ensembl Variant Effect Predictor (VEP).
# VEP is run in offline mode using the local cache and GRCh38 reference assembly.
# The output is generated as a tab-delimited annotation table.

set -euo pipefail


# ___ INPUT AND OUTPUT ___

INPUT_VCF="/path/to/input.vcf.gz"
OUT_TSV="/path/to/output/VEP_annotation.tsv"


# ___ TOOL ___

VEP=vep


# ___ RUN VEP ___
# Annotate variants using the local VEP cache.
# The output includes gene symbols, canonical transcripts, HGVS notation,
# transcript biotype, variant class, allele frequency information,
# and one selected annotation per variant (--pick).

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
    --pick

echo "Done: $OUT_TSV"