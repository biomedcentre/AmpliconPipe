#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly AMPLICON_REF=$1
readonly INPUT_BAM=$2 
readonly THREADS=$3
readonly PLOIDY=$4
readonly CONTAINER="$5"

readonly prefix=$(basename "${INPUT_BAM%%.sorted.dedup.bam}")
readonly output="/pipeline/output/${prefix}"

gatk HaplotypeCaller -R "${AMPLICON_REF}" -I "${INPUT_BAM}" -O "${output}".alt.ploidy.haplotype.caller.vcf.gz --native-pair-hmm-threads "${THREADS}" -ploidy "${PLOIDY}" --verbosity WARNING

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
if [ "${CONTAINER}" = 'docker' ]; then
    chown -Rc :"${GID:-0}" /pipeline/output
fi 
chmod -Rc g+w,o-rwx /pipeline/output