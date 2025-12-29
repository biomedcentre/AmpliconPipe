import pandas as pd 
import numpy as np
import logging

def read_vcf(link): 
    '''
    Reads vcf, extracts GT. Translates coords to hg38 coordinates so that resulting logs will contrain hg38 coordiantes  
    :param link: str, path to vcf file
    :return vcf: pd.Dataframe, vcf
    '''
    vcf = pd.read_csv(link, sep='\t', comment='#', header=None,
          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], index_col=False)  


    contig_name = vcf['CHROM'][0].split(':')[-1]
    if '-' in contig_name:
        contig_start = int(contig_name.split('-')[0])
    elif '+' in contig_name: 
        contig_start = int(contig_name.split('+')[0])
    else: 
        contig_start = int(contig_name)

     vcf = vcf.assign(CHROM=vcf['CHROM'].map(lambda x: x.split(':')[0]),
                    POS=(vcf['POS'] + contig_start - 1).astype(int),
                    GT=vcf['SAMPLE'].map(lambda x: x.split(':')[0]).str.replace('|', '/', regex=False))

    return vcf


def extract_genotypes(vcf):
    '''
    Transforms vcfs to dict of variant_ids POS_REF_ALT and their genotypes, vafs, ads, dps, quals 
    :param vcf: pd.DataFrame, read vcf 
    :return gt_vector: dict of lists 
    '''
    
    gt_vector = {}
    
    for idx, row in vcf.iterrows(): 
        if ',' in row['ALT']: # multiallelic sites are split into separate ids 

            for num, alt in enumerate(row['ALT'].split(',')):
                num = num+1
                variant_id = str(row['POS']) + '_' + row['REF'] + "_" + alt
                gt_relabled = row['GT'].split('/')
                gt_relabeled = [int(g == str(num)) for g in gt_relabled]
                gt_vector[variant_id] = '/'.join(gt_relabeled).strip('/')
        else:
            variant_id = str(row['POS']) + '_' + row['REF'] + "_" + row['ALT']
            gt_vector[variant_id] = gt
            
    return gt_vector


def genotype_split(genotype): 
    '''
    Splits genotype e.g 0/1 into int list 
    :param genotype: str, genotype 
    '''
    return [int(i) for i in genotype.split('/')]


def phase_variant(patient, mother, father, phased_variants, variant_name): 
    '''
    Phase variant with both parents 
    :param patient: str, patient genotype 
    :param mother: str, mother genotype 
    :param father: str, father genotype 
    :param phased_variants: dict, dict of lists of phased variants for trio, where keys are variant ids 
    :param variant_name: str, variant id 
    :return phased_variants: dict, updated dict 
    '''
    # transform genotypes 
    patient, mother, father = genotype_split(patient), genotype_split(mother), genotype_split(father)
    patient_state, mother_state, father_state = sum(patient), sum(mother), sum(father) 
    patient_phased, mother_phased, father_phased = {}, {}, {}

    # detect impossible states 
    min_patient_state = (mother_state // 2) + (father_state // 2)
    max_patient_state = (mother_state // 2) + (father_state // 2) + (mother_state % 2) + (father_state % 2)

    if patient_state > max_patient_state: 
        logging.warning(f'{variant_name} genotype in child is {patient}, parents are {father}, {mother}. Variant is de novo or one of the genotypes is incorrect. Additional checks are required')
        return phased_variants

    if patient_state < min_patient_state: 
        logging.warning(f'{variant_name} genotype in child is {patient}, parents are {father}, {mother}. One of the genotypes is incorrect. Additional checks are required')
        return phased_variants

    # detect unphasable or phasing not required variants 
    if_phased = True
    if (patient_state == mother_state) & (patient_state == father_state): 
        if_phased = False

    elif (patient_state == 0) & ((mother_state == 1) | (father_state == 1)):  # none of variants inherited by patient 
        patient_phased['M'], patient_phased['F'] = 0, 0 
        mother_phased['I'], father_phased['I'] = 0, 0 # none in inherited gaplotype 
        mother_phased['N'], father_phased['N'] = mother_state, father_state

    elif (patient_state == 1) & ((mother_state + father_state) == 1): 
        # patient has het, only mother or only father has het 
        patient_phased['M'], patient_phased['F'] = mother_state, father_state 
        mother_phased['I'], father_phased['I'] = mother_state, father_state 
        mother_phased['N'], father_phased['N'] = 0, 0 
        
    elif (patient_state == 1) & ((mother_state + father_state) > 1):
        # patient has het, mother or father has homo, other may have het
        patient_phased['M'], patient_phased['F'] = mother_state // 2, father_state // 2 # will be 1 where homo, 0 where parent has hetero or 0 
        mother_phased['I'], mother_phased['N'] = mother_state // 2, max(mother_state // 2, mother_state % 2)
        father_phased['I'], father_phased['N'] = father_state // 2, max(father_state // 2, father_state % 2)

     elif (patient_state == 2): 
        # patient has homo, parents have homo or het 
        patient_phased['M'], patient_phased['F'] = 1, 1 
        mother_phased['I'], mother_phased['N'] = 1, mother_state // 2
        father_phased['I'], father_phased['N'] = 1, father_state // 2

    # if phasing sucessful, update table of phased variants 
    if if_phased: 
        phased_variants[variant_name] = patient_phased['M'], patient_phased['F'], mother_phased['I'], mother_phased['N'], father_phased['I'], father_phased['N']

    return phased_variants
    
    
def phase_variant_duo(patient, parent, phased_variants, variant_name): 
    '''
    Phasing for duo with proxy as second parent. Proxy assumption is 1) 0/1 for every variant in proxy; 2) if patient has a variant and it is found in non-proxy parent, it assumed to be inherited from non proxy parent 
    :param patient: str, patient genotype 
    :param parent: str, patent genotype 
    :param phased_variants: dict, dict of lists of phased variants for trio, where keys are variant ids 
    :param variant_name: str, variant id 
    :return phased_variants: dict, updated dict 
    '''
    # transform genotypes 
    patient, parent, proxy = genotype_split(patient), genotype_split(parent), genotype_split('0/1')
    patient_state, parent_state, proxy_state = sum(patient), sum(parent), sum(proxy) 
    patient_phased, parent_phased, proxy_phased = {}, {}, {}

    # detect impossible states 
    min_patient_state = (parent_state // 2) # since min proxy is always 0 
    max_patient_state = (parent_state // 2) + (parent_state % 2) + 1 # max proxy is always 1 
    if (patient_state > max_patient_state) | (patient_state < min_patient_state): 
        logging.warning(f'{variant_name} genotype in child is {patient}, parents are {parent}, and proxy is 0/1. One of the genotypes is incorrect. Additional checks are required')
        return phased_variants

    # due to assumption there are no unphasable variants  

    if (patient_state == 0):  # none of variants inherited by patient, proxy is always het, homo parent is already filtered out  
        patient_phased['PA'], patient_phased['PR'] = 0, 0 
        parent_phased['I'], proxy_phased['I'] = 0, 0 # none in inherited gaplotype 
        parent_phased['N'], proxy_phased['N'] = parent_state, proxy_state

    elif (patient_state == 1) & (parent_state > 0): # parent has variant that patient has in het  
        # patient has het, only mother or only father has het 
        patient_phased['PA'], patient_phased['PR'] = 1, 0 
        parent_phased['I'], proxy_phased['I'] = 1, 0
        parent_phased['N'], proxy_phased['N'] = parent_state // 2, 1
        
    elif (patient_state == 1) & (parent_state == 0): # parent do not have variant from patient 
        patient_phased['PA'], patient_phased['PR'] = 0, 1 
        parent_phased['I'], proxy_phased['I'] = 0, 1
        parent_phased['N'], proxy_phased['N'] = 0, 0

    elif (patient_state == 2): # homozygous in patient 
        patient_phased['PA'], patient_phased['PR'] = 1, 1 
        parent_phased['I'], proxy_phased['I'] = 1, 1
        parent_phased['N'], proxy_phased['N'] = parent_state // 2, 0
        
    phased_variants[variant_name] = patient_phased['PA'], patient_phased['PR'], parent_phased['I'], parent_phased['N'], proxy_phased['I'], proxy_phased['N']

    return phased_variants 


def select_inhereted_h_from_paritialy_phased(patient, phased_mother, phased_father): 
    # take all variants already phased 
    # rephase 
    # find matches 


def phase_variants(): 

if __name__ == '__main__':  

    parser = argparse.ArgumentParser(description='Phases variants based on relatedness ')
    parser.add_argument('--input_folder', type=str, help='Input folder with vcfs', required=True)
    parser.add_argument('--vcf_postfix', type=str, help='Vcf postfix', default='.conflict.resolved.vcf.gz')
    parser.add_argument('--ped_file', type=str, help='Prepared ped file with relatedness', required=True)
    parser.add_argument('--output_prefix', type=str, help='Output_prefix', required=True)
    parser.add_argument('--gid', type=str, default='0')
    parser.add_argument('--container', type=str, default='singularity')
    args = parser.parse_args()

    ped_file = pd.read_csv(args.ped_file, sep='\t', header=None, index_col=None)

    # find all samples
    samples = []
    for line in ped_file.to_numpy(): 
        for s in line: 
            samples.append(s)
    samples = list(set(samples)) # add removal of "proxy" sample 

    # read data 
    vcfs, genotypes_extracted = {}, {}
    all_variants = []
    for s in samples: 
        vcfs[s] = read_vcf(os.path.join(input_folder, f'{s}{args.vcf_postfix}'))
        genotypes_extracted[s] = extract_genotypes(vcfs[s])
        all_variants = all_variants + list(genotypes_extracted[s].keys())

    all_variants = list(set(all_variants))

    phased_haplotypes = {s:pd.DataFrame([], columns=['hap1', 'hap2']) for s in samples}
    for idx, row in ped_file.iterrows(): 
        patient_id, mother_id, father_id = row[0], row[1], row[2] 
        if (mother_id != 'proxy') & (father_id != 'proxy'):
            # check if any haplotypes of mother and father are already phased 
            if (len(phased_haplotypes[mother_id]) > 0) | (len(phased_haplotypes[father_id]) > 0): 
                # write procedure for haplotype selection here 
                select_inhereted_h_from_paritialy_phased()
                # determine inherited for mother and father 
                # phase only variants in patient for variants in parents that are already phased 
                
            phased_variants = {}
            for variant_id in genotypes_extracted[patient_id]: 
                if variant_id not in all_variants: 
                    phased_variants = phase_variant(genotypes_extracted[patient_id][variant_id], genotypes_extracted[mother_id][variant_id], genotypes_extracted[father_id][variant_id], phased_variants, variant_id)

        else: 
            # proxy present 
            # check if any haplotypes of mother and father are already phased 
            if mother_id != 'proxy': 
                parent_id = mother_id
            else: 
                parent_id = father_id
            
            if len(phased_haplotypes[parent_id]) > 0: 
                # write procedure for haplotype selection here 
                select_inhereted_h_from_paritialy_phased()
                # determine inherited for mother and father 
                # phase only variants in patient for variants in parents that are already phased 
                
            phased_variants = {}
            for variant_id in genotypes_extracted[patient_id]: 
                if variant_id not in all_variants: 
                    phased_variants = phase_variant_duo(genotypes_extracted[patient_id][variant_id], genotypes_extracted[parent_id][variant_id], phased_variants, variant_name)