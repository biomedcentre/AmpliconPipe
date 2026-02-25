import os
import argparse
import logging
import sys
import subprocess 
import pandas as pd 
import numpy as np

components_location = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.join(components_location, 'resolve_conflicts'))

from resolve_conflicts import get_vcf_comments, genotype_vectors, write_vcf, restore_orig_coordinates, form_comments_for_final_vcf

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
                    POS=vcf['POS'] + contig_start - 1,
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
        group_id = str(row['POS']) + '_' + row['REF']
        
        if ',' in row['ALT']: # multiallelic sites are split into separate ids 
            for num, alt in enumerate(row['ALT'].split(',')):
                num = num+1
                variant_id = str(row['POS']) + '_' + row['REF'] + "_" + alt

                gt_relabled = row['GT'].split('/')
                gt_relabeled = [str(int(g == str(num))) for g in gt_relabled]
                gt_vector[variant_id] = [group_id, '/'.join(gt_relabeled).strip('/')]
        else:
            variant_id = str(row['POS']) + '_' + row['REF'] + "_" + row['ALT']
            gt_vector[variant_id] = [group_id, row['GT']]
    
    return pd.DataFrame(gt_vector, index=['group', 'gt']).T


def genotype_split(genotype): 
    '''
    Splits genotype e.g 0/1 into int list 
    :param genotype: str, genotype 
    '''
    return [int(i) for i in genotype.split('/')]


def return_genotype(variant_id, genotype_dict, is_deletion_source=False, is_patient_with_deletion=False): 
    '''
    Utility that returns genotype from dict if it is present, and 0/0 genotype if none 
    '''
    if variant_id in genotype_dict: 
        if (not is_patient_with_deletion) & (not is_deletion_source): 
            return genotype_dict[variant_id]
        elif is_patient_with_deletion: 
            if genotype_dict[variant_id] == '1/1':
                return '0/1' 
            else: 
                return genotype_dict[variant_id] # exclusion for 0/0 genotypes 
        elif is_deletion_source:
            return '0/0' # no variants considered inherited if deletion is none 
    else: 
        return '0/0'


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
        if (patient_state == 2) & ((mother_state > 0) | (father_state > 0)): # de novo homozygous with likely appearance
            patient_phased['M'], patient_phased['F'] = 1, 1
            mother_phased['I'], father_phased['I'] = max(mother_state // 2, mother_state % 2), max(father_state // 2, father_state % 2)
            mother_phased['N'], father_phased['N'] = mother_state % 2, father_state % 2     
            phased_variants[variant_name] = patient_phased['M'], patient_phased['F'], mother_phased['I'], mother_phased['N'], father_phased['I'], father_phased['N']
            return phased_variants 
        else:    
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
        
    elif (patient_state == 1) & ((mother_state == 2) | (father_state == 2)):
        # patient has het, mother or father has homo, other may have het
        patient_phased['M'], patient_phased['F'] = mother_state // 2, father_state // 2 # will be 1 where homo, 0 where parent has hetero or 0 
        mother_phased['I'], mother_phased['N'] = mother_state // 2, max(mother_state // 2, mother_state % 2)
        father_phased['I'], father_phased['N'] = father_state // 2, max(father_state // 2, father_state % 2)

    elif (patient_state == 2): 
        # patient has homo, parents have homo or het 
        patient_phased['M'], patient_phased['F'] = 1, 1 
        mother_phased['I'], mother_phased['N'] = 1, mother_state // 2
        father_phased['I'], father_phased['N'] = 1, father_state // 2
         
    else: 
        if_phased = False

    # if phasing sucessful, update table of phased variants 
    if if_phased: 
        phased_variants[variant_name] = patient_phased['M'], patient_phased['F'], mother_phased['I'], mother_phased['N'], father_phased['I'], father_phased['N']

    return phased_variants


def phase_variant_multiallelic(patient, mother, father, phased_variants, unphased_variant_name, phased_variant_name):
    '''
    If one of the alleles of multialleleic variant is unresolved, resolve it using phasing info of the second one 
    '''
    # transform genotypes 
    patient, mother, father = genotype_split(patient), genotype_split(mother), genotype_split(father)
    patient_state, mother_state, father_state = sum(patient), sum(mother), sum(father) 
    patient_phased, mother_phased, father_phased = {}, {}, {}

    # detect impossible states 
    min_patient_state = (mother_state // 2) + (father_state // 2)
    max_patient_state = (mother_state // 2) + (father_state // 2) + (mother_state % 2) + (father_state % 2)

    # assign heterozygous according to phased 
    already_phased = phased_variants[phased_variant_name] 
    if (sum(already_phased) == 0) | (sum(already_phased) == 6): # all homozygous, these samples wont be affected
        return phased_variants

    # assign using heterozygous variants
    if_phased = True
    if (patient_state == 1) & ((already_phased[1] + already_phased[0]) == 1): 
        patient_phased['M'], patient_phased['F'] = already_phased[1], already_phased[0]
    if (mother_state == 1) & ((already_phased[2] + already_phased[3]) == 1): 
        mother_phased['I'], mother_phased['N'] = already_phased[3], already_phased[2]
    if (father_state == 1) & ((already_phased[4] + already_phased[5]) == 1): 
        father_phased['I'], father_phased['N'] = already_phased[5], already_phased[4]

    # assign using homo variants in unphased variant
    if (patient_state != 1): 
        patient_phased['M'], patient_phased['F'] = patient_state // 2, patient_state // 2
    if mother_state != 1: 
        mother_phased['I'], mother_phased['N'] = mother_state // 2, mother_state // 2
    if father_state != 1: 
        father_phased['I'], father_phased['N'] = father_state // 2, father_phased // 2

    if (len(father_phased) > 0) & (len(mother_phased) > 0) & (len(patient_phased) > 0): # all phased, early exit 
        phased_variants[unphased_variant_name] = patient_phased['M'], patient_phased['F'], mother_phased['I'], mother_phased['N'], father_phased['I'], father_phased['N']
        return phased_variants

    if patient_state < min_patient_state: # variant absent in patient, most probable configuration can be suggested in some cases 
        if (mother_state == 2) & (len(father_phased) > 0):
            patient_phased['M'], patient_phased['F'] = patient_state // 2, father_phased['I']
        elif (father_state == 2) & (len(mother_phased) > 0):
            patient_phased['M'], patient_phased['F'] =  mother_phased['I'], patient_state // 2
        else: 
            phased_variants.pop(phased_variant_name)
            logging.warning(f'{variant_name} genotype in child is {patient}, parents are {father}, {mother}. One of the genotypes is incorrect. Additional checks are required')
            return phased_variants
    
    if patient_state > max_patient_state: # de novo variant appearance; 2/1 anynot1 anyone1 was already handled by previous logic. Other can't be solved
        phased_variants.pop(phased_variant_name)
        logging.warning(f'{variant_name} genotype in child is {patient}, parents are {father}, {mother}. Variant is de novo or one of the genotypes is incorrect. Additional checks are required')
        return phased_variants
                
    # assign those that do no have heterozygous, using recieved phasing, assume that there are no de novo variants 
    if (len(mother_phased) > 0) & (len(patient_phased) == 0): 
        patient_phased['M'], patient_phased['F'] = mother_phased['I'], patient_state - mother_phased['I']
    elif (len(father_phased) > 0) & (len(patient_phased) == 0): 
        patient_phased['M'], patient_phased['F'] = patient_state - father_phased['I'], father_phased['I'], 
 
    if len(patient_phased) > 0: 
        if len(mother_phased) == 0: 
            mother_phased['I'], mother_phased['N'] = patient_phased['M'], mother_state - patient_phased['M']
        if len(father_phased) == 0: 
            father_phased['I'], father_phased['N'] = patient_phased['M'], father_state - patient_phased['M']
    else: 
        if_phased = False
        
    if if_phased: 
        phased_variants[unphased_variant_name] = patient_phased['M'], patient_phased['F'], mother_phased['I'], mother_phased['N'], father_phased['I'], father_phased['N']

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


def process_inhereted_tab_to_hap_tab(ind_tab, block_id): 
    '''
    Function that forms haplotype tab from phased variants tab and adds proper block id (future PS tag) 
    :param ind_tab: pd.DataFrame, phased variants frame for individual 
    :param block_id: int, number of block = number of line in ped file
    :return hap_tab: pd.DataFrame, ind_tab renamed and with block added  
    '''
    hap_tab = []

    for idx, row in ind_tab.iterrows():
        if (row.sum() == 0) | (row.sum() == 2):
            block = -1
        else: 
            block = block_id
        hap_tab.append([row.iloc[0], row.iloc[1], block])

    hap_tab = pd.DataFrame(hap_tab, index=ind_tab.index, columns=['hap1', 'hap2', 'block'])
    
    return hap_tab 


def haplotype_append(hap_tab, ind_phased_haplotype, block_id): 
    '''
    Takes already phased haplotypes for individual and assesses if new phased variants can be used to extend already present blocks. Useful when more than one generation present/children from different parents are present 
    :param hap_tab: pd.DataFrame, ind_tab renamed and with block added  
    :param ind_phased_haplotype: pd.DataFrame, haplotypes with previously phased variants 
    :param block_id: int, number of block = number of line in ped file
    :return ind_phased_haplotype: pd.DataFrame, updated haplotypes 
    '''
    # check presense of any phased variants in haplotypes 
    if (ind_phased_haplotype['block'].fillna(-1) != -1).any(): 
        already_phased = ind_phased_haplotype[ind_phased_haplotype['block'].fillna(-1) != -1]
        for block, al_ph in already_phased.groupby('block'): 
            intersected_vars = hap_tab.index[hap_tab['block'] != -1].intersection(al_ph.index)
            if len(intersected_vars) > 0: # TO DO find proper contstant 
                if (hap_tab.loc[intersected_vars]['hap1'] == al_ph.loc[intersected_vars]['hap1']).all():  # TO DO: find proper constant for match
                    hap1, hap2, final_block = 'hap1', 'hap2', block
                    break
                elif (hap_tab.loc[intersected_vars]['hap2'] == al_ph.loc[intersected_vars]['hap1']).all(): 
                    hap1, hap2, final_block = 'hap2', 'hap1', block
                    break 
        else: 
            hap1, hap2, final_block = 'hap1', 'hap2', block_id 
    else: 
        hap1, hap2, final_block = 'hap1', 'hap2', block_id 

    for variant_id, row in hap_tab.iterrows():
        ind_phased_haplotype.loc[variant_id, hap1] = row['hap1']
        ind_phased_haplotype.loc[variant_id, hap2] = row['hap2']
        ind_phased_haplotype.loc[variant_id, 'block'] = final_block if row['hap1'] != row['hap2'] else -1 

    return ind_phased_haplotype


def assign_haplotypes_to_vcf(hap_tab, vcf): 
    '''
    Adds phasing info to vcf 
    :param hap_tab: pd.DataFrame, phased haplotypes
    :param vcf: pd.DataFrame, vcf read to pandas dataframe and processed 
    :return vcf: pd.DataFrame, vcf updated with phasing info 
    '''
    
    hap_tab = hap_tab[hap_tab['block'] != -1]

    if hap_tab.empty: 
        return vcf.drop(['GT'], axis=1)
    
    hap_tab['POS_REF'] = hap_tab.index.map(lambda x: '_'.join(x.split('_')[:-1]))
    vcf['POS_REF'] = vcf['POS'].astype(str) + '_' + vcf['REF']
    
    for pos_ref, tab in hap_tab.groupby('POS_REF'): 
        line_index = np.where(vcf['POS_REF'] == pos_ref)[0][0]
        if len(tab) == 1: 
            new_genotype = f"{tab.iloc[0]['hap1']}|{tab.iloc[0]['hap2']}"
        else:
            # take vcf row alt order, create mapping dict, create new genotype 
            map_dict = {alt:i+1 for i, alt in enumerate(vcf.loc[line_index]['ALT'].split(','))}
            tab['ALT'] = tab.index.map(lambda x: map_dict[x.split('_')[-1]])
            new_genotype = f"{(tab['ALT']*tab['hap1']).sum()}|{(tab['ALT']*tab['hap2']).sum()}"

        vcf.loc[line_index, 'SAMPLE'] = vcf.loc[line_index, 'SAMPLE'].replace(vcf.loc[line_index, 'GT'], new_genotype) + f":{tab['block'][0]}"
        vcf.loc[line_index, 'FORMAT'] = vcf.loc[line_index, 'FORMAT']  + ':PS'

        #vcf.loc[line_index] = vcf_in_pos_line

    return vcf.drop(['POS_REF', 'GT'], axis=1)


def form_comments_for_final_vcf(vcf_comments): 
    '''
    Forms new vcf header, combining headers of all used vcfs 
    :param vcf_comments: dict, comments from all 3 vcfs 
    :return comments: list, header of final vcf 
    '''
    comments = vcf_comments[:-1]
    vcf_col_names = vcf_comments[-1] 

    comments.append('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">\n') # add line about PS 
    comments.append(vcf_col_names)

    return comments


def ploidy_check(ids_list, ploidy):
    '''
    Function that checks ploidies in trio or duo and determines following actions 
    '''
    
    if all([ploidy[s] == 1 for s in ids_list]): # all have deletion, cis-trans phasing not applicable to any
        return 'continue'

    elif ploidy[ids_list[0]] != 1: # patient does not have deletion, one of parents may have deletion, phasing may proceeed as normal
        return 'normal phasing'

    elif (ploidy[ids_list[0]] == 1): # patient has deletion, need to understand, which copy is lost   
        return 'match copy strategy'

    else:
        return 'normal phasing'


def assess_mismatches(all_variants, patient_genotype, parent_genotype): 
    '''
    Function for patients with deletions. If patient has deletion, it determines from which parent copy came by finding parent 
    from whom the most variants in patient could be inherited.
    '''
    mismatches = 0
    for variant_id in all_variants: 
        patient_state = sum(genotype_split(return_genotype(variant_id, patient_genotype))) // 2 # 1 if hemizygous 
        parent_state = sum(genotype_split(return_genotype(variant_id, parent_genotype)))
        if patient_state > parent_state: 
            mismatches += 1
    return mismatches 


def determine_deletions_source(ids_list, genotypes_extracted, all_variants, ploidy):
    '''
    Function for patients with deletions. If patient has deletion, it determines from which parent copy came by finding parent 
    from whom the most variants in patient could be inherited. The parent from which copy was not inherited is deletion source/parent whose copy was lost
    '''
    if len(ids_list) == 2:  # deletion in duo, always assume it is inhereted 
        if ploidy[ids_list[1]] == 1: 
            deletion_source = ids_list[1] 
        else: 
            deletion_source = 'proxy'
    else:
        patient_id, mother_id, father_id = ids_list[0], ids_list[1], ids_list[2]
        mother_mismatch = assess_mismatches(all_variants, genotypes_extracted[patient_id]['gt'], genotypes_extracted[mother_id]['gt'])
        father_mismatch = assess_mismatches(all_variants, genotypes_extracted[patient_id]['gt'], genotypes_extracted[father_id]['gt'])
        
        if father_mismatch > mother_mismatch: 
            deletion_source = father_id
        elif father_mismatch < mother_mismatch:
            deletion_source = mother_id
        else: 
            if ploidy[father_id] == 1: 
                deletion_source = father_id
            elif ploidy[mother_id] == 1: 
                deletion_source = mother_id
            else: 
                logging.warning(f'Source of deleted copy in {patient_id} not identified. All parents have 2 copies and equal number of mistmatching variants. Randomly assigning father copy as lost copy')
                deletion_source = father_id

    return deletion_source 


def phase_trio_mode(patient_id, mother_id, father_id, genotypes_extracted, all_variants, deletion_source, is_patient_with_deletion): 
    '''
    Trio phasing 
    '''
    # TRIO PHASING 
    phased_variants = {}
    for variant_id in all_variants: 
        phased_variants = phase_variant(return_genotype(variant_id, genotypes_extracted[patient_id]['gt'], is_patient_with_deletion=is_patient_with_deletion),
                                            return_genotype(variant_id, genotypes_extracted[mother_id]['gt'], 
                                                            is_deletion_source=(mother_id == deletion_source)),
                                            return_genotype(variant_id, genotypes_extracted[father_id]['gt'],
                                                            is_deletion_source=(father_id == deletion_source)),
                                            phased_variants, variant_id)

    # check if any variants are multi alllelic 
    all_groups = pd.Series({v:'_'.join(v.split('_')[:-1]) for v in all_variants})
    multi_groups = all_groups.value_counts()[all_groups.value_counts() > 1]
    
    for group_id, count in multi_groups.items():
        v_in_group = all_groups[all_groups == group_id].index
        phased_alone = [v for v in v_in_group if v in phased_variants.keys()]
        if len(phased_alone) != 1: # no phased or both phased 
            continue 
        # if not, use multiallelic phasing to correct unphased variant
        phased, unphased = phased_alone[0], [v for v in v_in_group if v not in phased_variants.keys()][0]
        phased_variants = phase_variant_multiallelic(return_genotype(unphased, genotypes_extracted[patient_id]['gt'],
                                                                     is_patient_with_deletion=is_patient_with_deletion),
                                            return_genotype(unphased, genotypes_extracted[mother_id]['gt'], 
                                                            is_deletion_source=(mother_id == deletion_source)),
                                            return_genotype(unphased, genotypes_extracted[father_id]['gt'],
                                                            is_deletion_source=(father_id == deletion_source)), phased_variants, unphased, phased)   
        
    phased_variants = pd.DataFrame(phased_variants, index=['patient_1', 'patient_2', 'mother_1', 'mother_2', 'father_1', 'father_2']).T 
    
    return phased_variants 


def phase_duo_mode(patient_id, parent_id, genotypes_extracted, all_variants, deletion_source, is_patient_with_deletion):
    
    # phase variants with proxy    
    phased_variants = {}
    for variant_id in all_variants: 
            phased_variants = phase_variant_duo(return_genotype(variant_id, genotypes_extracted[patient_id]['gt'],
                                                                is_patient_with_deletion=is_patient_with_deletion),
                                                return_genotype(variant_id, genotypes_extracted[parent_id]['gt'],
                                                                is_deletion_source=(parent_id == deletion_source)), phased_variants, variant_id)

    # multiallelic is not required since all get resolved with proxy and inheritance from patent rule 
        
    phased_variants = pd.DataFrame(phased_variants, index=['patient_1', 'patient_2', 'parent_1', 'parent_2', 'proxy_1', 'proxy_2']).T 

    return phased_variants 
    

if __name__ == '__main__':  

    parser = argparse.ArgumentParser(description='Phases variants based on relatedness ')
    parser.add_argument('--input_folder', type=str, help='Input folder with vcfs', required=True)
    parser.add_argument('--vcf_postfix', type=str, help='Vcf postfix', default='.conflict.resolved.vcf.gz')
    parser.add_argument('--contamination_postfix', type=str, help='Vcf postfix', default='.contamination_ploidy_results_mqc.txt')
    parser.add_argument('--ped_file', type=str, help='Prepared ped file with relatedness', required=True)
    parser.add_argument('--family', type=str, help='Family id in ped file', required=True)
    parser.add_argument('--output_prefix', type=str, help='Output_prefix', required=True)
    parser.add_argument('--gid', type=str, default='0')
    parser.add_argument('--container', type=str, default='singularity')
    args = parser.parse_args()

    ped_file = pd.read_csv(args.ped_file, sep='\t', header=None, index_col=None)
    ped_file = ped_file[ped_file[0].astype(str) == args.family]

    # find all samples
    samples = []
    for line in ped_file.to_numpy(): 
        for s in line[1:-2]: 
            if s != 'proxy':
                samples.append(s)
    samples = list(set(samples)) 

    # read data and process vcfs into genotype dicts 
    vcfs, genotypes_extracted = {}, {}
    all_variants = []
    for s in samples: 
        vcfs[s] = read_vcf(os.path.join(args.input_folder, f'{s}{args.vcf_postfix}'))
        genotypes_extracted[s] = extract_genotypes(vcfs[s])
        all_variants = all_variants + list(genotypes_extracted[s].index)
        
    all_variants = list(set(all_variants))

    # read ploidy data
    ploidy = {}
    for s in samples: 
        ploidy[s] = pd.read_csv(os.path.join(args.input_folder, f'{s}{args.contamination_postfix}'), sep='\t')['ploidy'].iloc[0]
    ploidy['proxy'] = 3 # handle for processing proxy
    
    # modify ped file. replace parents that have 3 ploidy with proxy and drop patients with ploidy 3 
    drop_index = []
    for block_id, row in ped_file.iterrows():
        if (ploidy[row[1]] == 3) | ((ploidy[row[2]] == 3) & (ploidy[row[3]] == 3)): 
            drop_index.append(block_id)
            continue 
        if (ploidy[row[2]] == 3): 
            ped_file.loc[block_id, 2] = 'proxy'
            continue
        if (ploidy[row[3]] == 3): 
            ped_file.loc[block_id, 3] = 'proxy'

    ped_file =  ped_file.drop(drop_index)
    if ped_file.empty: 
        ped_file.to_csv(os.path.join(args.output_prefix, f'{args.family}.family.ped.txt'))
        sys.exit()
    
    # create empty tab of haplotypes 
    phased_haplotypes = {s:pd.DataFrame([], columns=['hap1', 'hap2', 'block']) for s in samples}
    
    # iterate over trios in ped file, phase each trio, merge haplotypes with haplotypes from previous trios if possible 
    for block_id, row in ped_file.iterrows(): 
        patient_id, mother_id, father_id = row[1], row[2], row[3] 
        
        if (mother_id != 'proxy') & (father_id != 'proxy'):
            placement_dict = {'patient': patient_id, 'mother': mother_id, 'father': father_id}
            # COPY NUMBER CHECKS 
            pl_ch = ploidy_check([patient_id, mother_id, father_id], ploidy)
        
            if pl_ch == 'continue': # if all have deletion, nothing to phase 
                continue
                
            if pl_ch == 'match copy strategy':
                # find deletion source and consider no variants inherited from deletion source 
                # phase considering variants from deletion source 0/0 and variants from patient 0/1 
                deletion_source = determine_deletions_source([patient_id, mother_id, father_id], genotypes_extracted, all_variants, ploidy)
                is_patient_with_deletion = True
            else: 
                deletion_source, is_patient_with_deletion = 'no', False
            
            phased_variants = phase_trio_mode(patient_id, mother_id, father_id, genotypes_extracted, all_variants, deletion_source, is_patient_with_deletion)
            
            # if deletion was present, remove deletion source and patient from haplotypes since they are not phasable anyway  
            if pl_ch == 'match copy strategy':
                drop_relative = 'mother' if deletion_source == mother_id else 'father'
                phased_variants = phased_variants.drop(['patient_1', 'patient_2', f'{drop_relative}_1', f'{drop_relative}_2'], axis=1)

        else: 
            # proxy present 
            parent_id = mother_id if mother_id != 'proxy' else father_id
            placement_dict = {'patient': patient_id, 'parent': parent_id}

            pl_ch = ploidy_check([patient_id, parent_id], ploidy)

            if pl_ch == 'continue': # if all have deletion, nothing to phase 
                continue
                
            if pl_ch == 'match copy strategy':
                # find deletion source and consider no variants inherited from deletion source 
                deletion_source = determine_deletions_source([patient_id, parent_id], genotypes_extracted, all_variants, ploidy)
                is_patient_with_deletion = True
            else: 
                deletion_source, is_patient_with_deletion = 'no', False

            phased_variants = phase_duo_mode(patient_id, parent_id, genotypes_extracted, all_variants, deletion_source, is_patient_with_deletion)
            
            # if deletion was present, remove deletion source and patient from haplotypes since they are not phasable anyway  
            if pl_ch == 'match copy strategy':
                drop_relative = 'parent' if deletion_source == parent_id else 'proxy'
                phased_variants = phased_variants.drop(['patient_1', 'patient_2', f'{drop_relative}_1', f'{drop_relative}_2'], axis=1)
        
        
        # if none phased, do not try to match haplotypes, just continue
        if phased_variants.empty: 
            continue

        for placement in phased_variants.columns.map(lambda x: x.split('_')[0]).unique():
            if placement != 'proxy':
                id_tab = process_inhereted_tab_to_hap_tab(phased_variants[[f'{placement}_1', f'{placement}_2']], block_id)
                # append to already present haplotypes, if possible, if no intersections, simply will be added to table as new block  
                phased_haplotypes[placement_dict[placement]] = haplotype_append(id_tab, phased_haplotypes[placement_dict[placement]], block_id)

    
    # modify vcfs with phased genotype and PS tag     
    for s in samples: 
        orig_link = os.path.join(args.input_folder, f'{s}{args.vcf_postfix}')
        vcf_comments = form_comments_for_final_vcf(get_vcf_comments(orig_link))
        final_vcf = assign_haplotypes_to_vcf(phased_haplotypes[s], vcfs[s])
        final_vcf = restore_orig_coordinates(final_vcf, orig_link)
        write_vcf(vcf_comments, final_vcf, f'{args.output_prefix}/{s}.family.phased.vcf')

    # temp output save 
    ped_file.to_csv(os.path.join(args.output_prefix, f'{args.family}.family.ped.txt'))

    if args.container == 'docker':
        subprocess.run(['chown', '-Rc', f':{args.gid}', args.output_prefix], capture_output=True, text=True, check=True)
    subprocess.run(['chmod', '-Rc', 'g+w,o-rwx', args.output_prefix], capture_output=True, text=True, check=True)