#!/bin/bash

#source "/pipeline/tools/header-settings.sh"

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ##
## Command-line arguments
readonly REF_FASTA="$1"
readonly INPUT_BAM="$2"
readonly REGIONS_BED="$3"

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ##
readonly OUTPUT_BASENAME="$(basename "${INPUT_BAM%%.*}")"
readonly VCF="${OUTPUT_BASENAME}.deepvariant.vcf"
readonly GVCF="${OUTPUT_BASENAME}.deepvariant.gvcf"

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ##
/opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref="${REF_FASTA}" \
  --reads="${INPUT_BAM}" \
  --regions="${REGIONS_BED}" \
  --output_vcf="/pipeline/output/${VCF}" \
  --num_shards="${THREADS:-23}" \
  --logging_dir="/pipeline/output/${OUTPUT_BASENAME}_deep_variant_logs" \
  --dry_run=false

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
chown -Rc :"${GID:-0}" /pipeline/output
chmod -Rc g+w,o-rwx /pipeline/output