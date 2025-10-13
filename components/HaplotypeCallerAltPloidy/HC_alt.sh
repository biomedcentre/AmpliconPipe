#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly AMPLICON_REF=$1
readonly BAM=$2 
readonly OUTPUT_PREFIX=$3
readonly THREADS=$4
readonly PLOIDY=$5

readonly prefix=$(basename "${INPUT_BAM%%.sorted.dedup.bam}")
readonly output="/pipeline/output/${prefix}"

gatk HaplotypeCaller -R "${AMPLICON_REF}" -I "${INPUT_BAM}" -O "${output}".haplotype.caller.altploidy.vcf.gz --native-pair-hmm-threads "${THREADS}" --verbosity WARNING -ploidy "${PLOIDY}"