import pandas as pd
import numpy as np
import argparse


def read_vcf(link, mode): 
    vcf = pd.read_csv(link, sep='\t', comment='#', header=None,
          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], index_col=False)  
    
    contig_start = int(vcf['CHROM'][0].split(':')[-1].split('+')[0]) # написать функции детекции сепараторов или паттерны 
    vcf = vcf.assign(CHROM=vcf['CHROM'].map(lambda x: x.split(':')[0]),
                    POS=vcf['POS'] + contig_start - 1,
                    GT=vcf['SAMPLE'].map(lambda x: x.split(':')[0]).str.replace('|', '/', regex=False))
    
    # add DP detection 
    if mode == 'dv': 
        vcf = vcf[vcf['FILTER'] == 'PASS']
        vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[1]))
    elif mode == 'hc': 
        vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[3]))
    elif mode == 'pileup': 
        vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[])) # узнать, где оно в пайлапе 
    
    return vcf


def get_vcf_comments(vcf_path):
    '''
    Extracts comments with header from vcf 
    :param vcf_path: str, link to vcf
    :return comments: list, all comments 
    '''
    comments = []
    
    if vcf_path.endswith('.gz'):
        with gzip.open(vcf_path, 'r') as f:
            a = f.readline().decode()
            while (a.startswith('#')):
                comments.append(a)
                a = f.readline().decode()
    else: 
        with open(vcf_path, 'r') as f:
            a = f.readline()
            while (a.startswith('#')):
                comments.append(a)
                a = f.readline()
    
    return comments


def write_vcf(comments, vcf_content, out_vcf):
    '''
    Writes vcf 
    :param comments: list, all comments 
    :param vcf_content: pd.DataFrame, vcf content
    :param out: str, link to out vcf
    '''
    
    with open(out_vcf , 'w') as f:
        for c in comments:
            f.write(c)
    vcf_content.to_csv(out_vcf , sep='\t', header=None, index=None, mode='a')
    
    
def genotype_vectors(vcf):
    
    gt_vector = {}
    
    for idx, row in vcf.iterrows():
        dp = row['DP']
        gt = gt.split('|')
        
        if ',' in row['ALT']:
            # function to check multiallelic indels like TTTG will be here 

            for num, alt in enumerate(row['ALT'].split(',')):
                num = num+1
                variant_id = str(row['POS']) + '_' + row['REF'] + "_" + alt
                
                vaf = int(row['AD'].split(',')[num])/row['DP']
                gt_val = sum([(g == str(num) for g in gt])  # this is required to work with ploidies more than 2 
                gt_vector[variant_id] = [gt_val, vaf, int(row['AD'].split(',')[num]),  dp, row['QUAL']]

        else:
            variant_id = str(row['POS']) + '_' + row['REF'] + "_" + row['ALT']
            gt_val = sum([(g == '1' for g in gt]) # get exact number of alleles to work with ploidy more than 2 
            vaf = int(row['AD'].split(',')[1])/row['DP']
            gt_vector[position] = [gt_val, vaf,  int(row['AD'].split(',')[num]), dp, row['QUAL']]

    return gt_vector


def hc_qual(hc_stats): 
    # checks quality of hc variant stats
    ht_gt, hc_vaf, hc_dp = hc_stats[0], hc_stats[1], hc_stats[3]
    if hc_dp >= 30: # check if HC stats can be used. probably add also QUAL check in future 
        if (hc_gt == 2 and hc_vaf > 0.85) or (hc_gt == 1 and hc_vaf >= 0.2 and hc_vaf < 0.85): 
            return True
    return False


def compare_positions_present_in_both_callers(resolved, log_buffer, gt_vectors_all):
    
    # variants called by both callers 
    both_variants_list = list(gt_vectors_all['hc'].keys() & gt_vectors_all['dv'].keys())
    for variant_id in both_variants_list: 
        if gt_vectors_all['hc'][variant_id][0] == gt_vectors_all['dv'][variant_id][0]: 
            resolved[variant_id] = 'hc'
        else: 
            # variant present but genotype differs 
            if variant_id in gt_vectors_all['pileup'].keys():
                pileup_gt = gt_vectors_all['pileup']
                for caller_type in ['hc', 'dv']: 
                    if pileup_gt == gt_vectors_all[caller_type][variant_id][0]:
                        log_buffer.append(f"[CONFLICT] Conflict at {variant_id}; {caller_type} genotype was supported by mpileup data")
                        resolved[variant_id] = caller_type
                        break 
                else: 
                     # if genotype mismatch here, write function to check if pileup correct
                    log_buffer.append(f"[WARNING] Conflict at {variant_id}; All genotypes were not supported by pileup data")
            else: 
                # currently pileup-call output does not have indels, this is why this code. Later provide indels (take them from mpileup or other thresholder) 
                if hc_qual(gt_vectors_all['hc'][variant_id])
                    current_resolution = 'hc' 
                else: 
                    current_resolution = 'dv' 
                # LATER ADD ADDITIONAL CHECK FOR DV QUAL 
                log_buffer.append(f"[CONFLICT] Conflict at {variant_id}; No bcftools call data. {current_resolution} genotype chosen due to HC good/bad qual")
                resolved[variant_id] = current_resolution
                
    return resolved, log_buffer


def check_positions_present_in_one_caller(resolved, log_buffer, gt_vectors_all, caller_type=caller_type)

    if caller_type == 'hc': 
        other = 'dv'
    else: 
        other = 'hc' 
        
    variants_to_check = set(gt_vectors_all[caller_type]) - set(gt_vectors_all[other])
    for variant_id in variants_to_check: 
        if caller_type == 'hc': # quality check, if qual okay, deem present
            if hc_qual(hc_stats): 
                resolved[variant_id] = 'hc'
                log_buffer.append(f"[CONFLICT] Conflict at {variant_id}; HC genotype chosen due to HC good/bad qual")
                continue
        # else go check pileups 
        if variant_id in gt_vectors_all['pileup'].keys(): 
            if gt_vectors_all[caller_type][variant_id][0] == gt_vectors_all['pileup'][variant_id][0]: # check if GT matches 
                resolved[variant_id] = caller_type
                log_buffer.append(f"[CONFLICT] Conflict at {variant_id}; {caller_type} genotype chosen due to pileup match")
            else: 
                # what to do if pileup genotype does not match; for MVP gonna trust in caller but later will change that 
                resolved[variant_id] = 'pileup' 
                log_buffer.append(f"[WARNING] Conflict at {variant_id}; Pileup genotype chosen; mismatch between callers")
        else:
            # need additional checks for indels here due to  indels not called now 
            log_buffer.append(f"[WARNING] Conflict at {variant_id}; absent in mpileup, dropping this variant")
        
    return resolved, log_buffer        


def compare_and_resolve(gt_vectors_all):
    resolved = dict()
    log_buffer = []

    # filter  vcfs for calling region earlier if required ADD EARLIER 

    # variants present in both caller 
    resolved, log_buffer = compare_positions_present_in_both_callers(resolved, log_buffer, gt_vectors_all)

    # variants in only one caller 
    for caller_type in ['hc', 'dv']:
        resolved, log_buffer = compare_positions_present_in_one_caller(resolved, log_buffer, gt_vectors_all, caller_type=caller_type)

    return resolved, log_buffer


def resolved_to_frame(resolved): 
    resolved_frame = []
    for k in resolved: 
        data = k.split('_')
        data.append(resolved_k)
        resolved_frame.append(data)
    resolved_frame = pd.DataFrame(resolved_frame, columns=resolved.keys(), index=['POS', 'REF', 'ALT', 'caller_type']).T

    return resolved_frame


def form_final_vcf(resolved, vcfs, genotype_vectors): 

    final_vcf = []
    for pos, all_tab in resolved.groupby('POS'): 
        # current chosen scheme for calls from different callers: add lines as multialleleic, later if needed, join with bcftools norm 
        # need to check DV formats though if merging is possible and needed 
        for caller_type, tab in tab.groupby('caller_type'):
            vcf_in_pos_line = vcfs[caller_type]['POS'] == pos].iloc[0]
            
            if len(vcf_in_pos_line['ALT'].split(',')) == len(tab): 
                new_genotype, quality_source = vcf_in_pos_line['GT'], '/'.join([caller_type] * ploidy)
            else: 
                # generate new genotype and replace old one; not removing alt alleles for now; might be removing if needed in future 
                new_genotype, quality_source, old_genotype = [], [], vcf_in_pos_line['GT'].split('/') 
                for num, alt in enumerate(vcf_in_pos_line['ALT'].split(',')): 
                    if (tab['ALT'] == alt).any(): 
                        gen_add, source = num+1, caller_type
                    else: 
                        gen_add, source = 0, 'pileup' 
                    for i in range(sum([g == num+1 for g in old_genotype])): 
                        new_genotype.append(gen_add)
                        quality_source.append(source)
                for i in range(sum([g == 0 for g in old_genotype])):  
                    new_genotype.append(0)
                    quality_source.append(caller_type)
                    
            vcf_in_pos_line = vcf_in_pos_line.assign(SAMPLE=vcf_in_pos_line['SAMPLE'].str.replace(vcf_in_pos_line['GT'], '/'.join(new_genotype)) + \
                                                     ':' + '/'.join(quality_source),
                                                     FORMAT=vcf_in_pos_line['FORMAT'] + ':SR')
            
            final_vcf.append(vcf_in_pos_line)    

    final_vcf = pd.concat(final_vcf, axis=1).drop(['GT', 'DP', 'AD'], axis=1).sort_values(by='POS')

    return final_vcf
            

if __name__ == '__main__':  

    parser = argparse.ArgumentParser(description='Compare variants between different callers for single sample')
    parser.add_argument('--hc_vcf', type=str, help='HaplotypeCaller vcf', required=True)
    parser.add_argument('--dv_vcf', type=str, help='DeepVariant vcf', required=True)
    parser.add_argument('--pileup', type=str, help='Bcftools pileup', required=True)
    parser.add_argument('--output_prefix', type=str, help='Output_prefix', required=True)
    args = parser.parse_args()
    print("PROGRAM ARGUMENTS: ", vars(args))

    vcf_links = {'hc': args.hc_vcf, 'dv': args.dv_vcf, 'pileup': args.pileup}
    gt_vectors, vcfs = {}, {}

    for key_name in vcf_links: 
        vcf = read_vcf(vcf_links[key_name], key_name)
        gt_vectors_all[key_name], vcfs[key_name] = genotype_vectors(vcf), vcf

    resolved, log_buffer = compare_and_resolve(gt_vectors_all)

    # write logs 
    with open(f'{args.output_prefix}.variant_conflict_resolution.log', 'w') as f: 
        for i in log_buffer:
            f.write(i + '\n')

    resolved = resolved_to_frame(resolved)
    final_vcf = form_final_vcf(resolved, vcfs, genotype_vectors)

    comments = get_vcf_comments(vcf_links['hc']) # MVP uses hc comments, that will be changed later to full updated comments 
    write_vcf(comments, final_vcf, f'{args.output_prefix}.conflict.resolved.vcf')