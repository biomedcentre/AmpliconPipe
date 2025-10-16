#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

## Command-line arguments
readonly REFERENCE="$1"

bwa-mem2/bwa-mem2 index "${REFERENCE}"

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
chown -Rc :"${GID:-0}" /pipeline/reference
chmod -Rc g+w,o-rwx /pipeline/reference