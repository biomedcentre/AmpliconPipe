#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

while getopts 'P:o:R:j:a:' flag
do
    case "${flag}" in
        P) PSEUDO_COORDS=${OPTARG};;
        o) PREFIX=${OPTARG};;
        R) FULL_REFERENCE=${OPTARG};;
        j) THREADS=${OPTARG};;
        a) AMPLICON_REF=${OPTARG};;
    esac
done

# generate pseudogene sequence 
samtools faidx "${FULL_REFERENCE}" "${PSEUDO_COORDS}" -o /pipeline/output/"${PREFIX}".temp.pseudo_seq.fasta
samtools faidx /pipeline/output/"${PREFIX}".temp.pseudo_seq.fasta

# USING NEAT 3.4 https://github.com/ncsa/NEAT/archive/refs/tags/3.4.zip
cd /pipeline/tools/NEAT-3.4
python "${PATH_TO_NEAT}"/gen_reads.py -r "${PREFIX}".temp.pseudo_seq.fasta -R 150 -o /pipeline/output/"${PREFIX}".temp.pseudo -c 100 --pe 150 30 -M 0 --rng 42 -E 0