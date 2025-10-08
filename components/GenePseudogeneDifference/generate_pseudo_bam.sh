#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

while getopts 'p:j:a:' flag
do
    case "${flag}" in
        p) GEN_FASTQ_PREFIX=${OPTARG};;
        j) THREADS=${OPTARG};;
        a) AMPLICON_REF=${OPTARG};;
    esac
done

# simulate variants from pseudogenic reads
bwa-mem2 mem -t "${THREADS}" "${AMPLICON_REF}" /pipeline/input/"${GEN_FASTQ_PREFIX}".temp.pseudo_read1.fq.gz /pipeline/input/"${GEN_FASTQ_PREFIX}".temp.pseudo_read2.fq.gz -U 6 |  # -U 4 -B 1 -O 3
samtools view -@ "${THREADS}" -S -b - |
samtools sort  -@ "${THREADS}" - > /pipeline/output/"${GEN_FASTQ_PREFIX}".pseudoalign.bam