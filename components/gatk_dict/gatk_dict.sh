#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly REFERENCE="$1"
readonly REF_DICT="$2"

gatk CreateSequenceDictionary --REFERENCE "${REFERENCE}" --OUTPUT "${REF_DICT}"

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
chown -Rc :"${GID:-0}" /pipeline/reference
chmod -Rc g+w,o-rwx /pipeline/reference