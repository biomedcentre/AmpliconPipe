#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset


while getopts 'F:o:R:j:N:a:' flag
do
    case "${flag}" in
        F) PSEUDO_COORDS=${OPTARG};;
        o) PREFIX=${OPTARG};;
        R) FULL_REFERENCE=${OPTARG};;
        j) THREADS=${OPTARG};;
        N) PATH_TO_NEAT=${OPTARG};;
        a) AMPLICON_REF=${OPTARG};;
    esac
done

# generate pseudogene sequence 
samtools faidx "${FULL_REFERENCE}" "${PSEUDO_COORDS}" -o "${PREFIX}".temp.pseudo_seq.fasta

# generation of simulated reads with NEAT
# TO DO: move to newer neat version or other read generator, ideally easily compartible with gatk
cd "${PATH_TO_NEAT}"
python "${PATH_TO_NEAT}"/gen_reads.py -r "${PREFIX}".temp.pseudo_seq.fasta -R 150 -o "${PREFIX}"/temp.pseudo -c 100 --pe 150 30 -M 0 

# simulate variants from pseudogenic reads
bwa-mem2 mem -t "${THREADS}" "${AMPLICON_REF}" "${PREFIX}"/temp.pseudo_read1.fq.gz "${PREFIX}"/temp.pseudo_read2.fq.gz -U 6 |  # -U 4 -B 1 -O 3
samtools view -@ "${THREADS}" -S -b - |
samtools sort  -@ "${THREADS}" - > "${PREFIX}"/pseudoalign.bam

gatk AddOrReplaceReadGroups \
                    -I "${PREFIX}"/pseudoalign.bam \
                    -O "${PREFIX}"/pseudoalign.replacerg.bam \
                    --RGSM "${PREFIX}" \
                    --RGLB lib1 \
                    --RGPL 'ILLUMINA' \
                    --RGPU test --VERBOSITY WARNING
                    
samtools index -@ "${THREADS}" "${PREFIX}"/pseudoalign.replacerg.bam
 
# remove intermediate files 
rm "${PREFIX}"/pseudoalign.bam \
   "${PREFIX}"/temp.pseudo_read1.fq.gz "${PREFIX}"/temp.pseudo_read2.fq.gz

# TO DO: write amplicon name extraction 
 
gatk HaplotypeCaller  \
         -R "${AMPLICON_REF}" \
         -I "${PREFIX}"/pseudoalign.replacerg.bam \
         -O "${PREFIX}"/AMPNAME_"${PSEUDO_COORDS}".difference.vcf \ 
         --native-pair-hmm-threads "${THREADS}" \
         --verbosity WARNING

rm "${PREFIX}"/pseudoalign.replacerg.bam "${PREFIX}"/pseudoalign.replacerg.bam.bai "${PREFIX}".temp.pseudo_seq.fasta