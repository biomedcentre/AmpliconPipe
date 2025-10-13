#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly AMPLICON_REF=$1
readonly INPUT_BAM=$2 
readonly OUTPUT_PREFIX=$3
readonly THREADS=$4
readonly PLOIDY=$5

readonly prefix=$(basename "${INPUT_BAM%%.sorted.dedup.bam}")
readonly output="/pipeline/output/${prefix}"

bcftools mpileup --threads "${THREADS}" -f "${AMPLICON_REF}" -L 20000 -a FORMAT/AD,FORMAT/DP --no-BAQ -d 3000 -Ou "${INPUT_BAM}" | \
bcftools call --ploidy "${PLOIDY}" -mv -Oz -o "${output}".naive.vcf.gz

#modifying VCF file (required to get rid of pool in header. if sm-ngp has this already fixed, it is not required) 
bcftools view -h "${output}".naive.vcf.gz > "${output}".header.txt
sed -i "s/pool/$prefix/" "${output}".header.txt
bcftools reheader -h "${output}".header.txt -o "${output}".bcftools.altploidy.vcf.gz "${output}".naive.vcf.gz

tabix "${output}".bcftools.altploidy.vcf.gz

rm "${output}".header.txt "${output}".naive.vcf.gz
