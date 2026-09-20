#!/usr/bin/env bash
# This script annotates a VCF file using SpliceAI.
# SpliceAI is run using the GRCh38 reference assembly.
# The output VCF is compressed and indexed after annotation.

set -euo pipefail


# ___ INPUT AND OUTPUT ___

INPUT_VCF="/path/to/input.vcf.gz"
OUTPUT_DIR="/path/to/output"
LOG_DIR="/path/to/logs"
REF="/path/to/reference/GRCh38.p14.genome.fa"


# ___ PARAMETERS ___

ASSEMBLY="grch38"
DISTANCE="500"
MASK="1"


# ___ TOOLS ___

SPLICEAI=spliceai
BGZIP=bgzip
TABIX=tabix


# ___ OUTPUT FILES ___

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

sample=$(basename "$INPUT_VCF" .vcf.gz)

OUT_VCF="${OUTPUT_DIR}/${sample}.spliceai.vcf"
STDOUT_LOG="${LOG_DIR}/${sample}.stdout.log"
STDERR_LOG="${LOG_DIR}/${sample}.stderr.log"


# ___ INPUT CHECKS ___

if [[ ! -f "$INPUT_VCF" ]]; then
    echo "ERROR: Input VCF not found" >&2
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference FASTA not found" >&2
    exit 1
fi


# ___ RUN SPLICEAI ___
# Annotate variants with SpliceAI using the GRCh38 reference genome.
# A maximum distance of 500 bp is considered around each variant.
# Masked scores are used for variant interpretation.

$SPLICEAI \
    -I "$INPUT_VCF" \
    -O "$OUT_VCF" \
    -R "$REF" \
    -A "$ASSEMBLY" \
    -D "$DISTANCE" \
    -M "$MASK" \
    >"$STDOUT_LOG" \
    2>"$STDERR_LOG"


# ___ COMPRESS AND INDEX OUTPUT ___

$BGZIP \
    -f \
    "$OUT_VCF"

$TABIX \
    -f \
    -p vcf \
    "${OUT_VCF}.gz"


echo "Done: ${OUT_VCF}.gz"