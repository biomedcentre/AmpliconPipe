#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

AMPLICON_REF=$1
BAM=$2 
OUTPUT_PREFIX=$3
THREADS=$4
SAMPLE_NAME=$5

bcftools mpileup --threads "${THREADS}" -f "${AMPLICON_REF}" -L 10000 -a FORMAT/AD,FORMAT/DP --no-BAQ -d 3000 -Ou "${BAM}" | \
bcftools call -mv -Oz -o "${OUTPUT_PREFIX}".naive.vcf.gz

#modifying VCF file (required to get rid of pool in header. if sm-ngp has this already fixed, it is not required) 
bcftools view -h "${OUTPUT_PREFIX}".naive.vcf.gz > "${OUTPUT_PREFIX}".header.txt
sed -i "s/pool/$SAMPLE_NAME/" "${OUTPUT_PREFIX}".header.txt
bcftools reheader -h "${OUTPUT_PREFIX}".header.txt -o "${OUTPUT_PREFIX}".bcftools.vcf.gz "${OUTPUT_PREFIX}".naive.vcf.gz

tabix "${OUTPUT_PREFIX}".bcftools.vcf.gz

rm "${OUTPUT_PREFIX}".header.txt "${OUTPUT_PREFIX}".naive.vcf.gz
