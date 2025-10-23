#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

## Command-line arguments
readonly REFERENCE="$1"
readonly CONTAINER="$2"

bwa-mem2/bwa-mem2 index "${REFERENCE}"

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
if [ "${CONTAINER}" = 'docker' ]; then
    chown -Rc :"${GID:-0}" /pipeline/output
fi 
chmod -Rc g+w,o-rwx /pipeline/reference