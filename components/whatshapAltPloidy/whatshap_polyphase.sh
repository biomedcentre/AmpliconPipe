#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly AMPLICON_REF=$1
readonly BAM=$2 
readonly FINAL_VCF=$3 # alt.ploidy.conflict.resolved.vcf
readonly OUTPUT_PREFIX=$4
readonly THREADS=$5
readonly PLOIDY=$6 

bgzip "${FINAL_VCF}"
tabix "${FINAL_VCF}".gz 
whatshap polyphase "${FINAL_VCF}".gz "${BAM}" -o "${OUTPUT_PREFIX}".phased.vcf --tag=PS -r "${AMPLICON_REF}" -p "${PLOIDY}" --threads "${THREADS}"
rm "${FINAL_VCF}".gz # remove unphased vcf 