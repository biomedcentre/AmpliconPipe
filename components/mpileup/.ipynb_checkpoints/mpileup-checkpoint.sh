#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

AMPLICON_REF=$1
BAM=$2 
FINAL_VCF=$3 
OUTPUT_PREFIX=$4
THREADS=$5
SAMPLE_NAME=$6

bcftools mpileup --threads  -f "${AMPLICON_REF}" -a FORMAT/AD,FORMAT/DP --no-BAQ -d 3000 -Ou "${BAM}" | \
bcftools call -mv -Oz -o "${OUTPUT_PREFIX}".naive.vcf.gz

#modifying VCF file (required to get rid of pool in header. if sm-ngp has this already fixed, it is not required) 
bcftools view -h "${file}" > "${OUTPUT_PREFIX}".header.txt
sed -i "s/pool/$SAMPLE_NAME/" "${OUTPUT_PREFIX}".header.txt

python "$current_dir/modify_vcf.py" "$file" temp.txt
cat header.txt temp.txt | bgzip > "$output_dir/${sample}.${type}.rehead.vcf.gz" && tabix "$output_dir/${sample}.${type}.rehead.vcf.gz"

rm "${OUTPUT_PREFIX}".header.txt "${OUTPUT_PREFIX}".temp.txt
