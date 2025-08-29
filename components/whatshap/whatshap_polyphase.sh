#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

AMPLICON_REF=$1
BAM=$2 
FINAL_VCF=$3 # alt.ploidy.conflict.resolved.vcf
OUTPUT_PREFIX=$4
THREADS=$5
PLOIDY=$6 

bgzip "${FINAL_VCF}"
tabix "${FINAL_VCF}".gz 
whatshap polyphase "${FINAL_VCF}".gz "${BAM}" -o "${OUTPUT_PREFIX}".phased.vcf --tag=PS -r "${AMPLICON_REF}" -p "${PLOIDY}" --threads "${THREADS}"
rm "${FINAL_VCF}".gz # remove unphased vcf 