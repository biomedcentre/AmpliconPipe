#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

AMPLICON_REF=$1
BAM=$2 
OUTPUT_PREFIX=$3
THREADS=$4

gatk HaplotypeCaller -R "${AMPLICON_REF}" -I "${BAM}" -O "${OUTPUT_PREFIX}".haplotype.caller.vcf.gz --native-pair-hmm-threads "${THREADS}" --verbosity WARNING