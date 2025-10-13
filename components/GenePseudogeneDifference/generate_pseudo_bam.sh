#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

while getopts 'P:j:a:' flag
do
    case "${flag}" in
        P) PSEUDO_COORDS=${OPTARG};;
        j) THREADS=${OPTARG};;
        a) AMPLICON_REF=${OPTARG};;
    esac
done

# simulate variants from pseudogenic reads
bwa-mem2/bwa-mem2 mem -t "${THREADS}" "${AMPLICON_REF}" /pipeline/input/"${PSEUDO_COORDS}".temp.pseudo_read1.fq.gz /pipeline/input/"${PSEUDO_COORDS}".temp.pseudo_read2.fq.gz -U 6 |  
samtools view -@ "${THREADS}" -S -b - |
samtools sort  -@ "${THREADS}" - > /pipeline/output/"${PSEUDO_COORDS}".pseudoalign.bam

rm /pipeline/output/"${PSEUDO_COORDS}".temp.pseudo_read1.fq.gz /pipeline/output/"${PSEUDO_COORDS}".temp.pseudo_read2.fq.gz 