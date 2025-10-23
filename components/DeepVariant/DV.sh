#!/bin/bash

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ##
## Command-line arguments
readonly REF_FASTA="$1"
readonly INPUT_BAM="$2"
readonly REGIONS_BED="$3"
readonly CONTAINER="$4"

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ##
readonly OUTPUT_BASENAME="$(basename "${INPUT_BAM%%.*}")"
readonly VCF="${OUTPUT_BASENAME}.deepvariant.vcf"
readonly GVCF="${OUTPUT_BASENAME}.deepvariant.gvcf"

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ##
DIR_TMP="/pipeline/output/${OUTPUT_BASENAME}.tmp.deepvariant"
mkdir -p "${DIR_TMP}"

## Attempt to redirect Bazel runfiles to a sample subdirectory
export TMPDIR="${DIR_TMP}"
export TEMP="${DIR_TMP}"
export TMP="${DIR_TMP}"

/opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref="${REF_FASTA}" \
  --reads="${INPUT_BAM}" \
  --regions="${REGIONS_BED}" \
  --output_vcf="/pipeline/output/${VCF}" \
  --num_shards="${THREADS:-23}" \
  --intermediate_results_dir="${DIR_TMP}" \
  --dry_run=false

rm -rdv --preserve-root "${DIR_TMP}"

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
if [ "${CONTAINER}" = 'docker' ]; then
    chown -Rc :"${GID:-0}" /pipeline/output
fi 
chmod -Rc g+w,o-rwx /pipeline/output