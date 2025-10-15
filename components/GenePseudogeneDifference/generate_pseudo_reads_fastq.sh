#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

while getopts 'P:R:j:' flag
do
    case "${flag}" in
        P) PSEUDO_COORDS=${OPTARG};;
        R) FULL_REFERENCE=${OPTARG};;
        j) THREADS=${OPTARG};;
    esac
done

# generate pseudogene sequence 
samtools faidx "${FULL_REFERENCE}" "${PSEUDO_COORDS}" -o /pipeline/output/"${PSEUDO_COORDS}".temp.pseudo_seq.fasta
samtools faidx /pipeline/output/"${PSEUDO_COORDS}".temp.pseudo_seq.fasta

# USING NEAT 3.4 https://github.com/ncsa/NEAT/archive/refs/tags/3.4.zip
cd /pipeline/tools/NEAT-3.4
python gen_reads.py -r /pipeline/output/"${PSEUDO_COORDS}".temp.pseudo_seq.fasta -R 150 -o /pipeline/output/"${PSEUDO_COORDS}".temp.pseudo -c 100 --pe 150 30 -M 0 --rng 42 -E 0

rm /pipeline/output/"${PSEUDO_COORDS}".temp.pseudo_seq.fasta /pipeline/output/"${PSEUDO_COORDS}".temp.pseudo_seq.fasta.fai

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
chown -Rc :"${GID:-0}" /pipeline/output
chmod -Rc g+w,o-rwx /pipeline/output