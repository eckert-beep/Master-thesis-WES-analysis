#!/usr/bin/env bash
# This script performs joint genotyping of multiple sample GVCF files.
# It checks and indexes individual GVCFs, combines them using GATK CombineGVCFs,
# and generates a jointly genotyped VCF using GATK GenotypeGVCFs.

# Safety settings:
# -e: exit on error
# -u: fail on unset variables
# -o pipefail: stop if any command in a pipeline fails
set -euo pipefail
shopt -s nullglob


# ___ INPUTS AND OUTPUTS ___
# Define input and output directories and reference files.

REF_GENOME=/path/to/reference/GRCh38.p14.genome.fa
GVCF_DIR=/path/to/gvcf
JOINT_DIR=/path/to/joint_genotyping

COMBINED_GVCF="${JOINT_DIR}/combined.g.vcf.gz"
JOINT_VCF="${JOINT_DIR}/joint_genotyped.vcf.gz"


# ___ THREADS AND TOOLS ___
# JAVA_MEM: maximum Java heap size for each GATK step

THREADS=64
GATK=gatk
TABIX=tabix

JAVA_MEM_COMBINE="-Xmx32g"     # 32 GB RAM for CombineGVCFs
JAVA_MEM_GENOTYPE="-Xmx24g"    # 24 GB RAM for GenotypeGVCFs

TMPDIR=/path/to/tmp             # temporary directory for intermediate files

mkdir -p "$JOINT_DIR" "$TMPDIR"


# ___ GVCF CHECK ___
# Check whether individual sample GVCF files are present.

if ! ls "${GVCF_DIR}"/*.g.vcf.gz >/dev/null 2>&1; then
    echo "No GVCF files found in ${GVCF_DIR}"
    exit 1
fi


# ___ GVCF INDEX CHECK ___
# Check whether all input GVCFs have a tabix index.
# Missing indexes are created automatically.

echo "Checking for missing .tbi indexes"

for GVCF in "${GVCF_DIR}"/*.g.vcf.gz; do

    if [[ ! -f "${GVCF}.tbi" ]]; then

        echo "Missing index: ${GVCF}.tbi"
        echo "Creating index..."

        $TABIX \
            -p vcf \
            "$GVCF"

        echo "Done."

    fi

done

echo "All GVCFs are indexed."


# ___ STEP 1: COMBINE GVCFs ___
# Combine all individual sample GVCFs into a single GVCF using GATK CombineGVCFs.
# Output: combined multi-sample GVCF

if [[ ! -f "$COMBINED_GVCF" ]]; then

    echo "Running CombineGVCFs on all GVCFs in: $GVCF_DIR"

    $GATK \
        --java-options "-Djava.io.tmpdir=${TMPDIR} ${JAVA_MEM_COMBINE}" \
        CombineGVCFs \
        -R "$REF_GENOME" \
        -O "$COMBINED_GVCF" \
        $(for GVCF in "${GVCF_DIR}"/*.g.vcf.gz; do
            echo -V "$GVCF"
        done)

else
    echo "Combined GVCF already exists, skipping CombineGVCFs: $(basename "$COMBINED_GVCF")"
fi


# ___ STEP 2: JOINT GENOTYPING ___
# Perform joint genotyping of the combined GVCF using GATK GenotypeGVCFs.
# Output: jointly genotyped multi-sample VCF

if [[ ! -f "$JOINT_VCF" ]]; then

    echo "Running GenotypeGVCFs on combined GVCF"

    $GATK \
        --java-options "-Djava.io.tmpdir=${TMPDIR} ${JAVA_MEM_GENOTYPE}" \
        GenotypeGVCFs \
        -R "$REF_GENOME" \
        -V "$COMBINED_GVCF" \
        -O "$JOINT_VCF"

else
    echo "Joint VCF already exists, skipping GenotypeGVCFs: $(basename "$JOINT_VCF")"
fi


# ___ STEP 3: INDEX JOINT VCF ___
# Create or overwrite the tabix index for the jointly genotyped VCF.

echo "Ensuring tabix index for joint VCF"

$TABIX \
    -f \
    -p vcf \
    "$JOINT_VCF"


echo
echo "Done"
echo "Combined GVCF: $COMBINED_GVCF"
echo "Joint VCF: $JOINT_VCF"