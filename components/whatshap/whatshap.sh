#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

AMPLICON_REF=$1
BAM=$2 
FINAL_VCF=$3 
OUTPUT_PREFIX=$4
THREADS=$5

whatshap phase "${FINAL_VCF}" "${BAM}" -o "${OUTPUT_PREFIX}".phased.vcf --tag=PS -r "${AMPLICON_REF}" --internal-downsampling 23 --threads "${THREADS}"