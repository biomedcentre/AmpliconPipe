#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

## Command-line arguments
readonly REF_PREFIX="$1"
readonly FASTQ_R1="$2"
readonly FASTQ_R2="$3"
readonly THREADS="$4"

readonly filename=$(basename "${FASTQ_R1%%.fastq.gz}")
sample_name="${filename%_L[0-9][0-9][0-9]_R[12]*}"
if [ "$sample_name" = "$filename" ]; then
    # If no lanes, try remove R group pattern only 
    sample_name="${filename%_R[12]*}"
fi

readonly OUTPUT_FASTP="/pipeline/output/${sample_name}.json"
readonly OUTPUT_BAM="/pipeline/output/${sample_name}.sorted.dedup.bam"

/pipeline/tools/fastp \
      -i "${FASTQ_R1}" \
      -I "${FASTQ_R2}" \
      --thread "${THREADS}" \
      --json="${OUTPUT_FASTP}" \
      --stdout \
    | bwa-mem2/bwa-mem2 mem \
      -p \
      -t "${THREADS}" \
      "${REF_PREFIX}" \
      -R "@RG\tID:sample\tPL:Illumina\tSM:${sample_name}" \
      /dev/stdin \
    | samtools view \
        -bh \
        --threads "${THREADS}" \
        - \
    | samtools fixmate \
        -m \
        --threads "${THREADS}" \
        - - \
    | samtools sort \
        --threads "${THREADS}" \
        -T '/pipeline/output/' \
        - \
    | samtools markdup \
        --threads "${THREADS}" \
        -T '/pipeline/output/' \
        -s \
        - \
        "${OUTPUT_BAM}"
        
samtools index \
    -@ "${THREADS}" \
    "${OUTPUT_BAM}"

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
chown -Rc :"${GID:-0}" /pipeline/output
chmod -Rc g+w,o-rwx /pipeline/output