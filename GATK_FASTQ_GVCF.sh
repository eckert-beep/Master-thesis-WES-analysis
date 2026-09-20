#!/usr/bin/env bash
# This script processes FASTQ files to generate GVCF files.
# It performs alignment, duplicate marking, base quality score recalibration,
# quality control, and variant calling.

# Safety settings:
# -e: exit on error
# -u: fail on unset variables
# -o pipefail: stop if any command in a pipeline fails
set -euo pipefail
shopt -s nullglob


# ___ INPUTS AND OUTPUTS ___
# Define input and output directories and reference files.

FASTQ_DIR=/path/to/fastq
REF_GENOME=/path/to/reference/GRCh38.p14.genome.fa
BAM_DIR=/path/to/bam
GVCF_DIR=/path/to/gvcf

# Known variant resources used for BQSR (Base Quality Score Recalibration)
DBSNP=/path/to/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz
MILLS=/path/to/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz


# ___ THREADS AND TOOLS ___
# THREADS: number of CPU cores per sample
# JAVA_MEM: maximum Java heap size for each GATK step

THREADS=32
BWA=bwa
SAMTOOLS=samtools
GATK=gatk
TABIX=tabix

JAVA_MEM_BQSR="-Xmx12g"     # 12 GB RAM for BQSR
JAVA_MEM_HC="-Xmx8g"        # 8 GB RAM for HaplotypeCaller

TMPDIR=/path/to/tmp          # temporary directory for intermediate files

mkdir -p "$BAM_DIR" "$GVCF_DIR" "$TMPDIR"


# ___ REFERENCE ___
echo "Ref checking indexes for $REF_GENOME"

[[ -f "${REF_GENOME}.fai" ]] || {
    echo "Ref create .fai"
    $SAMTOOLS faidx "$REF_GENOME"
}

DICT="${REF_GENOME%.*}.dict"

[[ -f "$DICT" ]] || {
    echo "Ref create .dict $DICT"
    $GATK CreateSequenceDictionary \
        -R "${REF_GENOME}" \
        -O "$DICT"
}

if [[ ! -f "${REF_GENOME}.bwt" && ! -f "${REF_GENOME}.0123" ]]; then
    echo "Ref create BWA index"
    $BWA index "$REF_GENOME"
fi


# ___ QC CHECK ___
qc () {
    local bam="$1"
    local tag="$2"

    [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" ]] || $SAMTOOLS index "$bam"

    $SAMTOOLS flagstat "$bam" > "${bam}.${tag}.flagstat.txt"
    $SAMTOOLS idxstats "$bam" > "${bam}.${tag}.idxstats.txt"
}


# ___ FASTQ SEARCH AND PAIR VALIDATION ___
# Loop through all sample folders in FASTQ_DIR.
# Each folder should contain paired-end FASTQ files:
#   *_1.fastq.gz  Read 1
#   *_2.fastq.gz  Read 2
#
# The script will:
#   - find all R1 files ending with "_1.fastq.gz"
#   - automatically detect matching R2 files
#   - skip incomplete pairs or interrupted transfers and issue a warning
#   - ignore FASTQC files

for PATDIR in "$FASTQ_DIR"/*/; do

    SAMPLE_ID=$(basename "$PATDIR")
    echo " Processing sample: $SAMPLE_ID "

    # Search for all *_1.fastq.gz files
    for R1 in "$PATDIR"/*_1.fastq.gz; do

        # Skip FASTQC files and incomplete transfers
        [[ "$R1" == *fastqc* || -e "$R1.part" ]] && continue

        base_name=$(basename "$R1")

        # Define matching R2 file by replacing "_1.fastq.gz" with "_2.fastq.gz"
        R2="${R1/%_1.fastq.gz/_2.fastq.gz}"

        if [[ ! -f "$R2" || -e "$R2.part" ]]; then
            echo "WARN: missing R2 for $(basename "$R1") --> skip"
            continue
        fi

        base_name=$(basename "$R1")             # NG-..._1.fastq.gz
        SAMPLE_BASE="${base_name%.fastq.gz}"    # NG-..._1
        OUT_BASE="${SAMPLE_BASE%_1}"            # NG-... (without _1/_2)


        # If both R1 and R2 exist, print the pair found
        echo " pair found:"
        echo "  R1: $(basename "$R1")"
        echo "  R2: $(basename "$R2")"


        # Read group:
        # SM: sample
        # PL: platform
        # LB: library
        # PU: platform unit
        RG="@RG\tID:${OUT_BASE}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:${OUT_BASE}\tPU:FLOWCELLX.L1"

        BAM_SORTED="${BAM_DIR}/${OUT_BASE}.sorted.bam"
        BAM_DEDUP="${BAM_DIR}/${OUT_BASE}.dedup.bam"
        METRICS="${BAM_DIR}/${OUT_BASE}.markdup.metrics.txt"
        RECAL_TABLE="${BAM_DIR}/${OUT_BASE}.recal.table"
        BAM_BQSR="${BAM_DIR}/${OUT_BASE}.recal.bam"
        GVCF="${GVCF_DIR}/${OUT_BASE}.g.vcf.gz"


        # ___ STEP 1: ALIGNMENT ___
        # Align paired-end reads (R1/R2) to the reference genome using BWA-MEM.
        # Output: sorted BAM

        if [[ ! -f "$BAM_SORTED" ]]; then

            echo " Aligning reads and sorting BAM "

            $BWA mem \
                -t "$THREADS" \
                -R "$RG" \
                "$REF_GENOME" \
                "$R1" \
                "$R2" \
                | $SAMTOOLS sort \
                    -@ "$THREADS" \
                    -o "$BAM_SORTED" -

            $SAMTOOLS index \
                -@ "$THREADS" \
                "$BAM_SORTED"

        else
            echo "Sorted BAM already exists, skipping alignment: $(basename "$BAM_SORTED")"
        fi

        qc "$BAM_SORTED" "sorted"


        # ___ STEP 2: MARK DUPLICATES ___
        # Mark duplicate reads using GATK MarkDuplicates.
        # Output: duplicate-marked BAM + metrics file

        if [[ ! -f "$BAM_DEDUP" ]]; then

            echo " Marking duplicates "

            $GATK \
                --java-options "-Djava.io.tmpdir=${TMPDIR} ${JAVA_MEM_BQSR}" \
                MarkDuplicates \
                -I "$BAM_SORTED" \
                -O "$BAM_DEDUP" \
                -M "$METRICS" \
                --CREATE_INDEX true \
                --VALIDATION_STRINGENCY LENIENT

        else
            echo "Deduplicated BAM already exists, skipping MarkDuplicates: $(basename "$BAM_DEDUP")"
        fi

        qc "$BAM_DEDUP" "dedup"


        # ___ STEP 3: BASE RECALIBRATION + APPLY BQSR ___
        # Recalibrate base quality scores to correct systematic errors.
        # Requires known variant resources (dbSNP and Mills).
        # Output: recalibrated BAM

        if [[ ! -f "$RECAL_TABLE" ]]; then

            echo " Base Quality Score Recalibration (BQSR) "

            $GATK \
                --java-options "-Djava.io.tmpdir=${TMPDIR} ${JAVA_MEM_BQSR}" \
                BaseRecalibrator \
                -I "$BAM_DEDUP" \
                -R "$REF_GENOME" \
                --known-sites "$DBSNP" \
                --known-sites "$MILLS" \
                -O "$RECAL_TABLE"

        else
            echo "Recalibration table already exists, skipping BQSR: $(basename "$RECAL_TABLE")"
        fi


        if [[ ! -f "$BAM_BQSR" ]]; then

            echo " Applying BQSR "

            $GATK \
                --java-options "-Djava.io.tmpdir=${TMPDIR} ${JAVA_MEM_BQSR}" \
                ApplyBQSR \
                -R "$REF_GENOME" \
                -I "$BAM_DEDUP" \
                --bqsr-recal-file "$RECAL_TABLE" \
                -O "$BAM_BQSR" \
                --create-output-bam-index true

        else
            echo "Recalibrated BAM already exists, skipping ApplyBQSR: $(basename "$BAM_BQSR")"
        fi

        qc "$BAM_BQSR" "recal"


        # ___ STEP 4: VARIANT CALLING ___
        # Call SNPs and indels per sample in GVCF mode using GATK HaplotypeCaller.
        # Output: GVCF file

        if [[ ! -f "$GVCF" ]]; then

            echo " Calling variants with HaplotypeCaller "

            $GATK \
                --java-options "-Djava.io.tmpdir=${TMPDIR} ${JAVA_MEM_HC}" \
                HaplotypeCaller \
                -R "$REF_GENOME" \
                -I "$BAM_BQSR" \
                -ERC GVCF \
                -O "$GVCF" \
                --native-pair-hmm-threads "$THREADS"

            echo " tabix index "

            $TABIX \
                -f \
                -p vcf \
                "$GVCF"

        else
            echo "GVCF already exists, skipping HaplotypeCaller: $(basename "$GVCF")"
        fi

        echo "Done: ${SAMPLE_ID} / ${OUT_BASE}"
        echo

    done
done