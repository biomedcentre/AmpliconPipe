#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly AMPLICON_REF=$1
readonly BAM=$2 
readonly FINAL_VCF=$3 # conflict.resolved.vcf
readonly CONTAINER="$4"

readonly prefix=$(basename "${BAM%%.sorted.dedup.bam}")

bgzip -f "${FINAL_VCF}" -o "${FINAL_VCF}".gz 
tabix "${FINAL_VCF}".gz 
whatshap phase "${FINAL_VCF}".gz "${BAM}" -o /pipeline/output/"${prefix}".phased.vcf --tag=PS -r "${AMPLICON_REF}" --internal-downsampling 23  

bgzip /pipeline/output/"${prefix}".phased.vcf
tabix /pipeline/output/"${prefix}".phased.vcf.gz

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
if [ "${CONTAINER}" = 'docker' ]; then
    chown -Rc :"${GID:-0}" /pipeline/output
fi 
chmod -Rc g+w,o-rwx /pipeline/output