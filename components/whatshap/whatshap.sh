#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

AMPLICON_REF=$1
BAM=$2 
FINAL_VCF=$3 # conflict.resolved.vcf
OUTPUT_PREFIX=$4

bgzip "${FINAL_VCF}"
tabix "${FINAL_VCF}".gz 
whatshap phase "${FINAL_VCF}".gz "${BAM}" -o "${OUTPUT_PREFIX}".phased.vcf --tag=PS -r "${AMPLICON_REF}" --internal-downsampling 23 
rm "${FINAL_VCF}".gz # remove unphased vcf 