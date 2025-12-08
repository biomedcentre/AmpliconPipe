import pandas as pd
import numpy as np
import argparse
import gzip 
import subprocess 
import os 
import sys 

components_location = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.join(components_location, 'resolve_conflicts'))

from resolve_conflicts import read_vcf, get_vcf_comments, genotype_vectors, write_vcf, resolved_to_frame, restore_orig_coordinates, form_comments_for_final_vcf


def qual(hc_stats): # for ploidy 3
    # checks quality of hc variant stats
    # TO DO: make threshold values constants 
    hc_gt, hc_vaf, hc_dp = hc_stats[0], hc_stats[1], hc_stats[3]
    if hc_dp >= 30: # check if HC stats can be used. probably add also QUAL check in future 
        if ((hc_gt == 1) and (hc_vaf > 0.1) and (hc_vaf < 0.495)) or ((hc_gt == 2) and (hc_vaf >= 0.495) and (hc_vaf < 0.85)) or ((hc_gt == 3) and (hc_vaf >= 0.85)): 
            return True
    return False


def form_final_vcf(resolved, vcfs, log_buffer): 
    '''
    Form contents of final vcf based on conflict resolution 
    :param resolved: pd.DataFrame of good variants
    :param vcfs: dict, dict of vcfs from different caller 
    :return final_vcf: pd.DataFrame, contents of final vcf 
    '''
    final_vcf = []
    for pos, tab in resolved.groupby('POS'): 
        vcf_in_pos_line = vcfs['hc'][vcfs['hc']['POS'] == pos].iloc[0]
        if len(vcf_in_pos_line['ALT'].split(',')) == len(tab):  # all genotypes passed qc 
            final_vcf.append(vcf_in_pos_line)    
        else: 
            # change to have only good allele genotype remaining?
            allele_num = len(vcf_in_pos_line['ALT'].split(','))
            log_buffer.append(f"[WARNING] {allele_num - len(tab)} alleles out of {allele_num} at {pos} were low quality; all variants in position will be removed from the final output")
            
    if len(final_vcf) > 0:
        final_vcf = pd.concat(final_vcf, axis=1).T.drop(['GT', 'DP', 'AD'], axis=1).sort_values(by='POS').reset_index(drop=True)
    else: 
        final_vcf = pd.DataFrame([], columns=vcfs['hc'].columns[:-3])
        log_buffer.append("[WARNING] No high-quality variants found, vcf is empty")

    return final_vcf, log_buffer

    
def check_quality(gt_vectors_all):
    '''
    Compares callers outputs and resolves conflicts
    :param gt_vectors_all: dict, of gt_vector frames for all callers 
    :param caller_type: str, caller type where variants are present, dv or hc 
    :return resolved: dict, variant_id:chosen caller_type. Variant ids that are found good with source caller 
    :return log_buffer: list, updated list of logged string about encountered conflicts
    '''
    resolved = dict()
    log_buffer = []

    for variant_id in gt_vectors_all['hc']: 
        if qual(gt_vectors_all['hc'][variant_id]): 
            resolved[variant_id] = 'hc'
        else: 
            log_buffer.append(f"[WARNING] Low HaplotypeCaller quality at {variant_id}; variant will be removed from final output")

    return resolved, log_buffer

    
if __name__ == '__main__':  

    parser = argparse.ArgumentParser(description='Compare variants between different callers for single sample')
    parser.add_argument('--hc_vcf', type=str, help='HaplotypeCaller vcf', required=True)
    parser.add_argument('--output_prefix', type=str, help='Output_prefix', required=True)
    parser.add_argument('--gid', type=str, default='0')
    parser.add_argument('--container', type=str, default='singularity')
    args = parser.parse_args()
    print("PROGRAM ARGUMENTS: ", vars(args))

    vcf_links = {'hc': args.hc_vcf}
    gt_vectors_all, vcfs, vcf_comments = {}, {}, {}

    for key_name in vcf_links: 
        vcf, vcf_comments[key_name] = read_vcf(vcf_links[key_name], key_name), get_vcf_comments(vcf_links[key_name])
        gt_vectors_all[key_name], vcfs[key_name] = genotype_vectors(vcf), vcf

    resolved, log_buffer = check_quality(gt_vectors_all)

    # write logs 
    with open(f'{args.output_prefix}.alt.ploidy.quality.log', 'w') as f: 
        for i in log_buffer:
            f.write(i + '\n')

    resolved = resolved_to_frame(resolved)
    final_vcf, log_buffer = form_final_vcf(resolved, vcfs, log_buffer)
    final_vcf = restore_orig_coordinates(final_vcf, vcf_links['hc'])
    comments =  form_comments_for_final_vcf(vcf_comments)
    
    write_vcf(comments, final_vcf, f'{args.output_prefix}.alt.ploidy.quality.checked.vcf')

    if args.container == 'docker':
        subprocess.run(['chown', '-Rc', f':{args.gid}', '/pipeline/output'], capture_output=True, text=True, check=True)
    subprocess.run(['chmod', '-Rc', 'g+w,o-rwx', '/pipeline/output'], capture_output=True, text=True, check=True)