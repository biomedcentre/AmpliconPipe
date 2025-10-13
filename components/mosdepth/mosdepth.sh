#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly INPUT_BAM="$1"
readonly REGIONS_BED="$2"
readonly THREADS="$3"

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ##
readonly prefix=$(basename "${INPUT_BAM%%.sorted.dedup.bam}")
readonly output="/pipeline/output/${prefix}"

/mosdepth \
  --threads "${THREADS}" \
  --by "${REGIONS_BED}" \
  --thresholds "1,10,20,50,100" \
  --no-per-base \
  /pipeline/output/"${prefix}" \
  "${INPUT_BAM}"

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
chown -Rc :"${GID:-0}" /pipeline/output
chmod -Rc g+w,o-rwx /pipeline/output