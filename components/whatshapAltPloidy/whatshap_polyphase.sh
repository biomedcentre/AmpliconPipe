#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly AMPLICON_REF=$1
readonly BAM=$2 
readonly FINAL_VCF=$3 # alt.ploidy.conflict.resolved.vcf
readonly THREADS=$4
readonly PLOIDY=$5 

bgzip "${FINAL_VCF}"
tabix "${FINAL_VCF}".gz 
whatshap polyphase "${FINAL_VCF}".gz "${BAM}" -o "${OUTPUT_PREFIX}".alt.ploidy.phased.vcf --tag=PS -r "${AMPLICON_REF}" -p "${PLOIDY}" --threads "${THREADS}"

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
chown -Rc :"${GID:-0}" /pipeline/output
chmod -Rc g+w,o-rwx /pipeline/output