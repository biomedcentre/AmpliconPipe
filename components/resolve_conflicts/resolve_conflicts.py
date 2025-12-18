import pandas as pd
import numpy as np
import argparse
import gzip 
import subprocess 

def read_vcf(link, mode): 
    '''
    Reads vcf, extracts GT, DP, AD for given caller type. Translates coords to hg38 coordinates so that resulting logs will contrain hg38 coordiantes  
    :param link: str, path to vcf file
    :param mode: str, caller type
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
    
    # TO DO: add calling region filter or add it to previous pipeline stages 
    
    # add DP detection 
    if mode == 'hc': 
        vcf = vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[1]),
                  DP=vcf['SAMPLE'].map(lambda x: x.split(':')[2]).astype(int))
    elif mode == 'dv': 
        vcf = vcf[vcf['FILTER'] == 'PASS']
        vcf = vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[3]),
                  DP=vcf['SAMPLE'].map(lambda x: x.split(':')[2]).astype(int))
    elif mode == 'pileup': 
        # TO DO modify multiallelic indels
        new_refalt = pd.DataFrame(list(vcf.apply(lambda x: modify_indel_bcftools(x['REF'], x['ALT'], x['INFO']), axis=1)), columns=['REF', 'ALT'])
        vcf = vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[3]),
                  DP=vcf['SAMPLE'].map(lambda x: x.split(':')[2]).astype(int),
                  REF=new_refalt['REF'],
                  ALT=new_refalt['ALT'],
                  ) 
    
    return vcf


def modify_indel_bcftools(ref, alt, info):
    '''
    Bcftools call returns indel REF and ALT in different format from HC and DV. To ensure REF/ALT match, transform records.
    :param ref: str, ref value
    :param alt: str, alt value
    :param info: str, alt value
    :return: str, str: new reference and new alt values 
    '''
    if info.startswith('INDEL'):
        if any([i.isupper() for i in alt]): # if it is insertion 
            inserted = np.where([i.isupper() for i in alt])[0]
            new_alt = alt[:inserted[0]] + alt[inserted[0]:inserted[-1]+1]
            new_ref = ref[:inserted[0]]
        else: 
            new_alt, new_ref = alt[:-1], ref[:-1]
    
        return new_ref.upper(), new_alt.upper() 
    else: 
        return ref, alt

    
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
    '''
    Transforms vcfs to dict of variant_ids POS_REF_ALT and their genotypes, vafs, ads, dps, quals 
    :param vcf: pd.DataFrame, read vcf 
    :return gt_vector: dict of lists 
    '''
    
    gt_vector = {}
    
    for idx, row in vcf.iterrows(): 
        dp = row['DP']
        gt = row['GT'].split('/')
        
        if ',' in row['ALT']: # multiallelic sites are split into separate ids 

            for num, alt in enumerate(row['ALT'].split(',')):
                num = num+1
                variant_id = str(row['POS']) + '_' + row['REF'] + "_" + alt
                vaf = int(row['AD'].split(',')[num])/row['DP']
                gt_val = sum([g == str(num) for g in gt])  # this is required to work with ploidies more than 2 
                gt_vector[variant_id] = [gt_val, vaf, int(row['AD'].split(',')[num]), dp, row['QUAL']]

        else:
            variant_id = str(row['POS']) + '_' + row['REF'] + "_" + row['ALT']
            gt_val = sum([g == '1' for g in gt]) # get exact number of alleles to work with ploidy more than 2 
            vaf = int(row['AD'].split(',')[1])/row['DP']
            gt_vector[variant_id] = [gt_val, vaf, int(row['AD'].split(',')[1]), dp, row['QUAL']]
    
    return gt_vector


def hc_qual(hc_stats): # correct thresholds for ploidy 2 only, will need to change that 
    '''
    Check HC variant correspondence to set quality thresholds 
    :param hc_stats: list, list of stats for variant_id from gt_vector
    :return: bool, True if good, and False if bad
    '''
    
    # checks quality of hc variant stats
    # TO DO: make threshold values constants 
    hc_gt, hc_vaf, hc_dp = hc_stats[0], hc_stats[1], hc_stats[3]
    if hc_dp >= 30: # check if HC stats can be used. probably add also QUAL check in future 
        if (hc_gt == 2 and hc_vaf > 0.85) or (hc_gt == 1 and hc_vaf >= 0.2 and hc_vaf < 0.85): 
            return True
    return False


def compare_positions_present_in_both_callers(resolved, log_buffer, gt_vectors_all):
    '''
    Compares genotypes of variant ids present in both (HC&DV) callers. If the same, takes HC version, if different, resolves using pileup vcf 
    :param resolved: dict, variant_id:chosen caller_type
    :param log_buffer: list, list of logged string about encountered conflicts
    :param gt_vectors_all: dict, of gt_vector frames for all callers 
    :return resolved: dict, updated variant_id:chosen caller_type
    :return log_buffer: list, updated list of logged string about encountered conflicts
    '''
    # variants called by both callers 
    both_variants_list = list(gt_vectors_all['hc'].keys() & gt_vectors_all['dv'].keys())
    for variant_id in both_variants_list: 
        if gt_vectors_all['hc'][variant_id][0] == gt_vectors_all['dv'][variant_id][0]: # genotypes match, everything okay
            resolved[variant_id] = 'hc'
        else: 
            # variant present but genotype differs 
            if variant_id in gt_vectors_all['pileup'].keys():
                pileup_gt = gt_vectors_all['pileup']
                for caller_type in ['hc', 'dv']: 
                    if pileup_gt == gt_vectors_all[caller_type][variant_id][0]: # found matching genotype in bcftools 
                        log_buffer.append(['[CONFLICT]', variant_id, 'Genotype corrected', caller_type.upper(), 'Genotype different between callers. Genotype supported by bcftools data chosen'])
                        resolved[variant_id] = caller_type
                        break 
                else: 
                    resolved[variant_id] = 'pileup'  # LATER: write function to check if pileup correct
                    log_buffer.append(['[WARNING]', variant_id, 'Genotype corrected', 'pileup', 'All genotypes were not supported by bcftools data, choosing bcftools genotype'])
            else: 
                # bcftools calls indels, but may miss some, this is why this code 
                if hc_qual(gt_vectors_all['hc'][variant_id]):
                    current_resolution = 'hc' 
                else: 
                    # LATER ADD ADDITIONAL CHECK FOR DV QUAL 
                    current_resolution = 'dv' 
                log_buffer.append(['[CONFLICT]', variant_id, 'Genotype corrected', current_resolution.upper(), 'Genotype chosen with HC quality assesement'])
                resolved[variant_id] = current_resolution
                
    return resolved, log_buffer


def check_positions_present_in_one_caller(resolved, log_buffer, gt_vectors_all, caller_type='hc'):
    '''
    Checks validity of variants that are present in given caller type, but absent in other caller type
    :param resolved: dict, variant_id:chosen caller_type
    :param log_buffer: list, list of logged string about encountered conflicts
    :param gt_vectors_all: dict, of gt_vector frames for all callers 
    :param caller_type: str, caller type where variants are present, dv or hc 
    :return resolved: dict, updated variant_id:chosen caller_type
    :return log_buffer: list, updated list of logged string about encountered conflicts
    '''
    if caller_type == 'hc': 
        other = 'dv'
    else: 
        other = 'hc' 
        
    variants_to_check = set(gt_vectors_all[caller_type]) - set(gt_vectors_all[other])
    
    for variant_id in variants_to_check: 
        if caller_type == 'hc': # quality check, if qual okay, deem present
            if hc_qual(gt_vectors_all['hc'][variant_id]): 
                resolved[variant_id] = 'hc'
                log_buffer.append(['[CONFLICT]', variant_id, 'Accepted', 'HC', 'Variant mismatch between callers. HC version chosen due to HC good qual'])
                continue
        # else go check pileups; this may change if DV qual thresholds appear 
        if variant_id in gt_vectors_all['pileup'].keys(): 
            if gt_vectors_all[caller_type][variant_id][0] == gt_vectors_all['pileup'][variant_id][0]: # check if GT matches with pileup
                resolved[variant_id] = caller_type
                log_buffer.append(['[CONFLICT]', variant_id, 'Accepted', caller_type.upper(), f'Variant mismatch between callers. Variant present in bcftools, caller accepted.'])
            else: 
                # TO DO: what to do if pileup genotype does not match; for MVP gonna trust in pileup but later probably will change that 
                resolved[variant_id] = 'pileup' 
                log_buffer.append(['[WARNING]', variant_id, 'Accepted with different genotype', 'pileup', f'Variant mismatch between callers. Variant present in {caller_type}. Variant present in bcftools. Bcftools genotype chosen'])
        else:
            # TO DO: some indels may be still absent in bcftools, see if resolvable
            log_buffer.append(['[WARNING]', variant_id, 'Discarded', 'pileup', f'Variant mismatch between callers. Variant present in {caller_type.upper()}. Variant absent in bcftools. Discarding this variant'])
        
    return resolved, log_buffer        


def compare_and_resolve(gt_vectors_all):
    '''
    Compares callers outputs and resolves conflicts
    :param gt_vectors_all: dict, of gt_vector frames for all callers 
    :param caller_type: str, caller type where variants are present, dv or hc 
    :return resolved: dict, variant_id:chosen caller_type. Variant ids that are found good with source caller 
    :return log_buffer: list, updated list of logged string about encountered conflicts
    '''
    resolved = dict()
    log_buffer = []

    # variants present in both caller 
    resolved, log_buffer = compare_positions_present_in_both_callers(resolved, log_buffer, gt_vectors_all)

    # variants in only one caller 
    for caller_type in ['hc', 'dv']:
        resolved, log_buffer = check_positions_present_in_one_caller(resolved, log_buffer, gt_vectors_all, caller_type=caller_type)

    return resolved, log_buffer


def resolved_to_frame(resolved): 
    '''
    :param resolved: dict, variant_id:chosen caller_type. Variant ids that are found good with source caller 
    :return resolved_frame: pd.DataFrame, resolved transformed to frame with split variant id  
    '''
    
    resolved_frame = []
    for k in resolved: 
        data = k.split('_')
        data.append(resolved[k])
        resolved_frame.append(data)
    resolved_frame = pd.DataFrame(resolved_frame, index=resolved.keys(), columns=['POS', 'REF', 'ALT', 'caller_type'])
    resolved_frame = resolved_frame.assign(POS=resolved_frame['POS'].astype(int))
    
    return resolved_frame


def form_final_vcf(resolved, vcfs): 
    '''
    Form contents of final vcf based on conflict resolution 
    :param resolved: pd.DataFrame of good variants
    :param vcfs: dict, dict of vcfs from different caller 
    :return final_vcf: pd.DataFrame, contents of final vcf 
    '''
    final_vcf = []
    for pos, all_tab in resolved.groupby('POS'): 
        # current chosen scheme for calls from different callers: add lines as multialleleic, later if needed, join with bcftools norm 
        # need to check DV formats though if merging is possible and needed 
        for caller_type, tab in all_tab.groupby('caller_type'):
            vcf_in_pos_line = vcfs[caller_type][vcfs[caller_type]['POS'] == pos].iloc[0]

            if len(vcf_in_pos_line['ALT'].split(',')) == len(tab): 
                # if all alts classified as from this caller, no changes in line are required 
                new_genotype, quality_source = vcf_in_pos_line['GT'], '/'.join([caller_type] * len(vcf_in_pos_line['GT'].split('/')))
            else: 
                # generate new genotype and replace old one; not removing alt alleles for now; might be removing if needed in future 
                new_genotype, quality_source, old_genotype = [], [], vcf_in_pos_line['GT'].split('/') 
                for num, alt in enumerate(vcf_in_pos_line['ALT'].split(',')): 
                    if (tab['ALT'] == alt).any(): 
                        gen_add, source = num+1, caller_type
                    else: 
                        gen_add, source = 0, 'pileup'           
                    for i in range(sum([int(g) == num+1 for g in old_genotype])): 
                        new_genotype.append(str(gen_add))
                        quality_source.append(source)
                for i in range(sum([g == 0 for g in old_genotype])):  
                    new_genotype.append(0)
                    quality_source.append(caller_type)
                new_genotype, quality_source = '/'.join(new_genotype), '/'.join(quality_source)
                
            vcf_in_pos_line['SAMPLE'] = vcf_in_pos_line['SAMPLE'].replace(vcf_in_pos_line['GT'], new_genotype) 
            
            # TO DO: find place to log quality source in vcf if needed 
            #+ ':' + quality_source 
            #vcf_in_pos_line['FORMAT'] = vcf_in_pos_line['FORMAT'] + ':SR'
            
            final_vcf.append(vcf_in_pos_line)    

    final_vcf = pd.concat(final_vcf, axis=1).T.drop(['GT', 'DP', 'AD'], axis=1).sort_values(by='POS').reset_index(drop=True)

    return final_vcf


def restore_orig_coordinates(final_vcf, orig_link):
    '''
    Originaly vcf had coordinates in amplicon ref. To use with bam file, transform to this coords back
    :param final_vcf: pd.DataFrame, contents of final vcf 
    :param orig_link: str, link to one of the original vcfs 
    :return final_vcf: pd.DataFrame, contents of final vcf 
    '''
    # if vcf had coordinates translated 

    orig_vcf = pd.read_csv(orig_link, sep='\t', comment='#', header=None,
          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], index_col=False)  

    contig_name = orig_vcf['CHROM'][0].split(':')[-1]
    if '-' in contig_name:
        contig_start = int(contig_name.split('-')[0])
    elif '+' in contig_name: 
        contig_start = int(contig_name.split('+')[0])
    else: 
        contig_start = int(contig_name)

    final_vcf = final_vcf.assign(CHROM=orig_vcf['CHROM'].iloc[0],
                          POS=(final_vcf['POS'].astype(int) - contig_start + 1))
    
    return final_vcf


def form_comments_for_final_vcf(vcf_comments): 
    '''
    Forms new vcf header, combining headers of all used vcfs 
    :param vcf_comments: dict, comments from all 3 vcfs 
    :return comments: list, header of final vcf 
    '''
    comments = vcf_comments['hc'][:-1]
    vcf_col_names = vcf_comments['hc'][-1] 

    for caller_type in ['dv', 'pileup']: 
        if caller_type in vcf_comments.keys():
            for line in vcf_comments[caller_type][:-1]: 
                line_start = line.split(',')[0]
                if not any([i.startswith(line_start) for i in comments]):
                    comments.append(line)
    comments.append(vcf_col_names)

    return comments 


if __name__ == '__main__':  

    parser = argparse.ArgumentParser(description='Compare variants between different callers for single sample')
    parser.add_argument('--hc_vcf', type=str, help='HaplotypeCaller vcf', required=True)
    parser.add_argument('--dv_vcf', type=str, help='DeepVariant vcf', required=True)
    parser.add_argument('--pileup', type=str, help='Bcftools pileup', required=True)
    parser.add_argument('--output_prefix', type=str, help='Output_prefix', required=True)
    parser.add_argument('--gid', type=str, default='0')
    parser.add_argument('--container', type=str, default='singularity')
    args = parser.parse_args()
    print("PROGRAM ARGUMENTS: ", vars(args))

    vcf_links = {'hc': args.hc_vcf, 'dv': args.dv_vcf, 'pileup': args.pileup}
    gt_vectors_all, vcfs, vcf_comments = {}, {}, {}

    for key_name in vcf_links: 
        vcf, vcf_comments[key_name] = read_vcf(vcf_links[key_name], key_name), get_vcf_comments(vcf_links[key_name])
        gt_vectors_all[key_name], vcfs[key_name] = genotype_vectors(vcf), vcf

    resolved, log_buffer = compare_and_resolve(gt_vectors_all)

    # write logs 
    log_buffer = pd.DataFrame(log_buffer, columns=['type', 'variant_id', 'verdict', 'correct_caller', 'details'])
    log_buffer.to_csv(f'{args.output_prefix}.variant_conflict_resolution.log', sep='\t', index=None)

    resolved = resolved_to_frame(resolved)
    final_vcf = form_final_vcf(resolved, vcfs)
    final_vcf = restore_orig_coordinates(final_vcf, vcf_links['hc'])
    comments =  form_comments_for_final_vcf(vcf_comments)
    
    write_vcf(comments, final_vcf, f'{args.output_prefix}.conflict.resolved.vcf')

    if args.container == 'docker':
        subprocess.run(['chown', '-Rc', f':{args.gid}', f'{args.output_prefix}.conflict.resolved.vcf', f'{args.output_prefix}.variant_conflict_resolution.log'], capture_output=True, text=True, check=True)
    subprocess.run(['chmod', '-Rc', 'g+w,o-rwx', f'{args.output_prefix}.conflict.resolved.vcf', f'{args.output_prefix}.variant_conflict_resolution.log'], capture_output=True, text=True, check=True)