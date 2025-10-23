#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly AMPLICON_REF=$1
readonly INPUT_BAM=$2 
readonly THREADS=$3
readonly CONTAINER="$4"

readonly prefix=$(basename "${INPUT_BAM%%.sorted.dedup.bam}")
readonly output="/pipeline/output/${prefix}"

bcftools mpileup --threads "${THREADS}" -f "${AMPLICON_REF}" -L 20000 -a FORMAT/AD,FORMAT/DP --no-BAQ -d 3000 -Ou "${INPUT_BAM}" | \
bcftools call -mv -Oz -o "${output}".bcftools.vcf.gz

tabix "${output}".bcftools.vcf.gz

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
if [ "${CONTAINER}" = 'docker' ]; then
    chown -Rc :"${GID:-0}" /pipeline/output
fi 
chmod -Rc g+w,o-rwx /pipeline/output