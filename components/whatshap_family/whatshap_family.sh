#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

readonly FAMILY_NAME=$1
readonly REHEADER=${2:-'missing'}
readonly REHEADER_BACK=${3:-'missing'}
readonly FAMILY_PHASING_FOLDER=$4
readonly BAM_LINKS=$5 # this will be a string of space separated files 
readonly VCF_LINKS=$6 # this will be a string of space separated files 
readonly CONTAINER=$7


tmp_folder=/pipeline/output/FAMILY_"${FAMILY_NAME}"_TMP
mkdir -p "${tmp_folder}"

# ==== REHEAD VCFS FOR MERGE IF NEEDED ==== 

if [ "${REHEADER}" != 'missing' ]; then 
    vcf_links=''
    vcfarray=($VCF_LINKS)
    for vcf in "${vcfarray[@]}"
        do
            # reheader vcfs so they can be merged 
            prefix=$(basename "${vcf%%.conflict.resolved.vcf.gz}")
            "${tmp_vcf}" = "${tmp_folder}"/"${prefix}".temp.modified.vcf.gz
            
            bcftools annotate --rename-chrs "${REHEADER}" -o "${tmp_vcf}" "${vcf}"
            tabix "${tmp_vcf}"
            vcf_links= "${vcf_links} ${tmp_vcf}"
        done
else 
    vcf_links="${VCF_LINKS}"
fi 

# ==== MERGE ALL VCF FOR WHATSHAP INPUT ==== 
    
if [ "${REHEADER_BACK}" = 'missing' ]; then
    bcftools merge --missing-to-ref -o "${tmp_folder}"/"${FAMILY_NAME}".merged.vcf "${vcfs}"
else 
    bcftools merge --missing-to-ref -o "${tmp_folder}"/"${FAMILY_NAME}".merged.v1.vcf "${vcfs}"
    bcftools annotate --rename-chrs "${REHEADER_BACK}" -o ${tmp_folder}/"${FAMILY_NAME}".merged.vcf "${tmp_folder}"/"${FAMILY_NAME}".merged.v1.vcf
fi 

# ==== WHATSHAP ==== 
output="${FAMILY_PHASING_FOLDER}"/FAMILY_"${FAMILY_NAME}"_phased.vcf
    
whatshap phase "${tmp_folder}"/"${FAMILY_NAME}".merged.vcf \
                            "${bam_links}" \
                            -o "${output}" \ 
                            --tag=PS -r "${REFERENCE}" --recombrate 0.01 \
                            --ped "${tmp_folder}"/"${FAMILY_NAME}".ped.txt

bgzip -f "${output}"
tabix "${output}".gz

# ==== TECH AND CLEAN UP ==== 
rm -rd "${tmp_folder}"

## These commands ensure that the output is not saved with root:root ownership.
## GID is group id env variable 
if [ "${CONTAINER}" = 'docker' ]; then
    chown -Rc :"${GID:-0}" /pipeline/output
fi 
chmod -Rc g+w,o-rwx /pipeline/output