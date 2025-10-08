#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

while getopts 'p:P:j:a:' flag
do
    case "${flag}" in
        p) GEN_FASTQ_PREFIX=${OPTARG};;
        P) PSEUDO_COORDS=${OPTARG};;
        j) THREADS=${OPTARG};;
        a) AMPLICON_REF=${OPTARG};;
    esac
done

gatk AddOrReplaceReadGroups \
                    -I /pipeline/input/"${GEN_FASTQ_PREFIX}".pseudoalign.bam \
                    -O /pipeline/output/"${GEN_FASTQ_PREFIX}".pseudoalign.replacerg.bam \
                    --RGSM "${GEN_FASTQ_PREFIX}" \
                    --RGLB lib1 \
                    --RGPL 'ILLUMINA' \
                    --RGPU test --VERBOSITY WARNING
                    
samtools index -@ "${THREADS}" /pipeline/output/"${GEN_FASTQ_PREFIX}".pseudoalign.replacerg.bam
 
AMPNAME=$(basename "${AMPLICON_REF}")
AMPNAME="${AMPNAME%.*}" 

gatk HaplotypeCaller \
         -R "${AMPLICON_REF}" \
         -I /pipeline/output/"${GEN_FASTQ_PREFIX}".pseudoalign.replacerg.bam \
         -O /pipeline/output/"${AMPNAME}".vs."${PSEUDO_COORDS}".difference.vcf \
         --native-pair-hmm-threads "${THREADS}" \
         --verbosity WARNING

rm /pipeline/output/"${GEN_FASTQ_PREFIX}".pseudoalign.replacerg.bam*