#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly ORIGINAL_NAME=$1 
readonly OUTPUT=$2 
readonly CONTAINER=$3

NEW_NAME=$(echo "${ORIGINAL_NAME}" | sed 's/:/_/g')

echo -e "${ORIGINAL_NAME} ${NEW_NAME}\n" > "${OUTPUT}"/reheader.txt
echo -e "${NEW_NAME} ${ORIGINAL_NAME}\n" > "${OUTPUT}"/reheader_back.txt

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
if [ "${CONTAINER}" = 'docker' ]; then
    chown -Rc :"${GID:-0}" /pipeline/output
fi 
chmod -Rc g+w,o-rwx /pipeline/output

