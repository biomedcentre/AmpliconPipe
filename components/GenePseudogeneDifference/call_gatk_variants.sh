#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

while getopts 'p:P:j:a:' flag
do
    case "${flag}" in
        P) PSEUDO_COORDS=${OPTARG};;
        j) THREADS=${OPTARG};;
        a) AMPLICON_REF=${OPTARG};;
    esac
done

gatk AddOrReplaceReadGroups \
                    -I /pipeline/input/"${PSEUDO_COORDS}".pseudoalign.bam \
                    -O /pipeline/output/"${PSEUDO_COORDS}".pseudoalign.replacerg.bam \
                    --RGSM "${PSEUDO_COORDS}" \
                    --RGLB lib1 \
                    --RGPL 'ILLUMINA' \
                    --RGPU test --VERBOSITY WARNING
                    
samtools index -@ "${THREADS}" /pipeline/output/"${PSEUDO_COORDS}".pseudoalign.replacerg.bam
 
AMPNAME=$(basename "${AMPLICON_REF}")
AMPNAME="${AMPNAME%.*}" 

gatk HaplotypeCaller \
         -R "${AMPLICON_REF}" \
         -I /pipeline/output/"${PSEUDO_COORDS}".pseudoalign.replacerg.bam \
         -O /pipeline/output/"${AMPNAME}".vs."${PSEUDO_COORDS}".difference.vcf \
         --native-pair-hmm-threads "${THREADS}" \
         --verbosity WARNING

rm /pipeline/output/"${PSEUDO_COORDS}".pseudoalign.replacerg.bam /pipeline/output/"${PSEUDO_COORDS}".pseudoalign.replacerg.bam.bai  

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
chown -Rc :"${GID:-0}" /pipeline/output
chmod -Rc g+w,o-rwx /pipeline/output