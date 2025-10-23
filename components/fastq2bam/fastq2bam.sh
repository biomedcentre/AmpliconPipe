#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

## Command-line arguments
readonly REF_PREFIX="$1"
readonly FASTQ_R1="$2"
readonly FASTQ_R2="$3"
readonly THREADS="$4"
readonly CONTAINER="$5"
readonly sample_name="$6" #snakemake sample wildcard

readonly OUTPUT_FASTP="/pipeline/output/${sample_name}.json"
readonly OUTPUT_BAM="/pipeline/output/${sample_name}.sorted.dedup.bam"

/pipeline/tools/fastp \
      -i "${FASTQ_R1}" \
      -I "${FASTQ_R2}" \
      --thread "${THREADS}" \
      --json="${OUTPUT_FASTP}" \
      --stdout \
    | /pipeline/tools/bwa-mem2/bwa-mem2 mem \
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

## These commands ensure that the output is not saved with root:root ownership in docker 
## GID is group id env variable 
if [ "${CONTAINER}" = 'docker' ]; then
    chown -Rc :"${GID:-0}" /pipeline/output
fi 
chmod -Rc g+w,o-rwx /pipeline/output