#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly AMPLICON_REF=$1
readonly BAM=$2 
readonly FINAL_VCF=$3 # alt.ploidy.conflict.resolved.vcf
readonly THREADS=$4
readonly PLOIDY=$5 
readonly CONTAINER="$6"

readonly prefix=$(basename "${BAM%%.sorted.dedup.bam}")

bgzip "${FINAL_VCF}" -o "${FINAL_VCF}".gz
tabix "${FINAL_VCF}".gz 

whatshap polyphase "${FINAL_VCF}".gz "${BAM}" -o /pipeline/output/"${prefix}".alt.ploidy.phased.vcf --tag=PS -r "${AMPLICON_REF}" -p "${PLOIDY}" --threads "${THREADS}" && \
bgzip /pipeline/output/"${prefix}".alt.ploidy.phased.vcf && \
tabix /pipeline/output/"${prefix}".alt.ploidy.phased.vcf.gz

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
if [ "${CONTAINER}" = 'docker' ]; then
    chown -Rc :"${GID:-0}" /pipeline/output
fi 
chmod -Rc g+w,o-rwx /pipeline/output