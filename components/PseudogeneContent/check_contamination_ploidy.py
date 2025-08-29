import argparse 
import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
from math import sqrt
from sklearn.mixture import GaussianMixture


# TO DO: extract constants from all functions and fix them here 

def read_vcf(link): 
    '''
    Reads vcf, extracts GT, DP, AD for given caller type. Translates coords to hg38 coordinates so that resulting logs will contrain hg38 coordiantes  
    :param link: str, path to vcf file
    :param mode: str, caller type
    :return vcf: pd.Dataframe, vcf
    '''
    vcf = pd.read_csv(link, sep='\t', comment='#', header=None,
          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], index_col=False)  
    
    contig_start = int(vcf['CHROM'][0].split(':')[-1].split('+')[0]) # TO DO: write more universal pattern function able to work with different contig namings
    vcf = vcf.assign(CHROM=vcf['CHROM'].map(lambda x: x.split(':')[0]),
                    POS=vcf['POS'] + contig_start - 1,
                    GT=vcf['SAMPLE'].map(lambda x: x.split(':')[0]).str.replace('|', '/', regex=False))

    # TO DO: add calling region filter or add it to previous pipeline stages 
    
    vcf = vcf.assign(AD=vcf['SAMPLE'].map(lambda x: x.split(':')[1]),
                  DP=vcf['SAMPLE'].map(lambda x: x.split(':')[2]).astype(int))
    
    return vcf


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

    gt_vector = pd.DataFrame(gt_vector, index=['GT', 'VAF', 'AD', 'DP', 'QUAL']).T
    
    return gt_vector
    

def choose_best_gm(baf_input, max_components=4): 
    '''
    Chooses best gaussian mixture for given baf input 
    :param baf_input: pd.Series, vector of bafs for sample
    :param max_ploidy: int, maximal ploidy allowed in fit 
    :return gm: sklearn.GaussianMixture, fitted model 
    '''
    baf_input = np.array(baf_input).reshape(-1, 1)
    
    # chose best gaussian
    bic = {} 
    for i in range(1, max_components+1): 
        gm = GaussianMixture(n_components=i, random_state=0).fit(baf_input) 
        bic[i] = gm.bic(baf_input)
        
    # choose best components 
    best_n = pd.Series(bic).idxmin()
    gm = GaussianMixture(n_components=best_n, random_state=0).fit(baf_input) 
    
    return best_n, gm 


def vafs_for_ploidy(ploidy):
    '''
    Generate vector of expected vaf peak means for given ploidy
    :param ploidy: int, ploidy
    :return vafs: list, list of floats, expected vafs 
    '''
    step = 1/ploidy
    vafs = []
    for i in range(1, ploidy+1): 
        vafs.append(step*i)

    return vafs


def plot_vaf_vector(vaf_vector, prefix, fit_type): 
    '''
    Plots current vaf vector
    :param vaf_vector: pd.Series, vector of vafs for variants
    :param prefix: str, path to output
    :param fit_type: str, description of current set of variants 
    '''
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    sns.distplot(vaf_vector, ax=ax)
    ax.set_title(f'{fit_type} VAF')
    ax.set_xlim(-0.05, 1.05)
    plt.rcParams['figure.dpi'] = 150
    plt.savefig(f'{prefix}.{fit_type}.VAF.png')


def concat_close_baf_peaks(good_peaks_stats, pred_labels, vaf_vector, allowed_std, max_ploidy):
    '''
    Gaussian mixture may create extra peaks that actually belong to same expected VAF peak. To correct that, close peaks will be concated.
    Closeness is calculated based on intersections between peaks (based on std calculation). 
    Extra threshold to avoid concating peaks that may intersect due to larger stds but are in different expected VAF peaks is introdused 
    :param good_peaks_stats: pd.DataFrame, mean and std of peaks
    :param pred_labels: pd.Series, predicted peaks for variants of vaf vector
    :param vaf_vector: pd.Series, vector of vafs for variants
    :param allowed_std: float, num of stds that will be used to calculate peak area 
    :param max_ploidy: int, maximal ploidy allowed in fit 
    :return new_good_peaks_stats: pd.DataFrame, updated peak stats 
    '''
    intersection_thres = 0.3
    peak_dist_thres = 1/(max_ploidy + 1)
    
    good_peaks_stats.index = good_peaks_stats['name']  
    peak_groups, current_group, current_stats = [], [good_peaks_stats['name'].iloc[0]], good_peaks_stats.iloc[0]
    
    for idx, row in good_peaks_stats.iloc[1:].iterrows():
        # extract peak intervals for groups 
        interval = (current_stats['mean'] - allowed_std*current_stats['std'], current_stats['mean'] + allowed_std*current_stats['std'])
        interval_row = (row['mean'] - allowed_std*row['std'], row['mean'] + allowed_std*row['std'])
        intersection_percentage = (interval_row[1] - interval[0])/(allowed_std*row['std']*2)

        # peaks that are too close get placed in one group 
        if (intersection_percentage > intersection_thres) & (abs(row['mean'] - current_stats['mean']) < peak_dist_thres): 
            current_group.append(idx)   
            # update stats
            current_group_samples = pred_labels.isin(current_group)
            current_stats = {'mean': vaf_vector[current_group_samples].mean(), 'std': vaf_vector[current_group_samples].std()}
        else:  # if next peak is too far save previous group and start anew 
            peak_groups.append(current_group)
            current_group, current_stats = [idx], good_peaks_stats.loc[idx]
    else: # last group in cycle  
        peak_groups.append(current_group)

    # form new peak stats 
    new_good_peaks_stats = []
    for group in peak_groups: 
        if len(group) == 1: # one peak in group - use its earlier stats 
            new_good_peaks_stats.append(list(good_peaks_stats.loc[group[0]]))
        else: 
            current_group_samples = pred_labels.isin(group)
            new_good_peaks_stats.append([''.join([str(g) for g in group]), vaf_vector[current_group_samples].mean(), vaf_vector[current_group_samples].std()])

    new_good_peaks_stats = pd.DataFrame(new_good_peaks_stats, columns=['name', 'mean', 'std'])
    
    return new_good_peaks_stats
            

def estimate_peak_fits(vaf_vector, log_buffer, prefix, max_ploidy=3, fit_type='All variants'):
    '''
    Fits GM for VAF vector and checks list of ploidies to see if any fits. 
    :param vaf_vector: pd.Series, vector of vafs for variants
    :param log_buffer: list, list of logged string about ploidy fits
    :param prefix: str, path to output
    :param max_ploidy: int, max ploidy that will be tested. Default is 3 for better peak resolution
    :param fit_type: str, description of current set of variants 
    :return: bool, if there was enough data for fit
    :return ploidy: bool or int, fit result
    :return log_buffer: list, updated log_buffer 
    '''
    minimal_peak_size = 5
    allowed_std = 3 

    plot_vaf_vector(vaf_vector, prefix, fit_type) # plotting at this point for now

    # choose best GM and filter before merge 
    best_n, gm = choose_best_gm(vaf_vector)
    pred_labels = pd.Series(gm.predict(np.array(vaf_vector).reshape(-1, 1)), index=vaf_vector.index)
    good_peaks = pred_labels.value_counts()[pred_labels.value_counts() >= minimal_peak_size].index
    
    if len(good_peaks) == 0: # no peaks with enough variants 
        log_buffer.append(f'[WARNING] {fit_type} fit: each peak in the best fit contains less than 5 variants. Contamination and ploidy will not be estimated')
        return False, False, log_buffer 
    
    good_peaks_stats = []
    for i in good_peaks: 
        mu, std = gm.means_[:,0][i], sqrt(gm.covariances_[:,0][:,0][i])
        good_peaks_stats.append([i, mu, std]) 

    good_peaks_stats = pd.DataFrame(good_peaks_stats, columns=['name', 'mean', 'std']).sort_values(by='mean', ascending=False)
    good_peaks_stats = concat_close_baf_peaks(good_peaks_stats, pred_labels, vaf_vector, allowed_std, max_ploidy).sort_values(by='mean', ascending=False)
 
    ploidy = False 
    for p in range(1, max_ploidy+1): # iterating through ploidies to see if any fits 
        expected_vafs = vafs_for_ploidy(p)
        
        if len(expected_vafs) < len(good_peaks_stats):
            log_buffer.append(f'[INFO] {fit_type} fit: Number of VAF peaks more than expected for ploidy {p}. Ploidy rejected')
            continue 
            
        peak_match = {}
        for idx, row in good_peaks_stats.iterrows(): # looking for matches for peaks of fit
            for match in expected_vafs[::-1]: # starting from larger vafs as stats table also starts from larger stats 
                if (match > (row['mean'] - allowed_std*row['std'])) & (match < (row['mean'] + allowed_std*row['std'])):  
                    peak_match[row['name']] = match 
                    expected_vafs = [i for i in expected_vafs if i < match]
                    break # match found, skip iteration over expected vafs
            else: 
                log_buffer.append(f'[INFO] {fit_type} fit: No match found for peak with mean {row["mean"]} and std ({row["std"]}) in VAF peaks expected for ploidy {p}. Ploidy rejected')
                break # further iteration not required, one mismatching peaks is enough
                
        if len(peak_match) == len(good_peaks_stats):
            log_buffer.append(f'[INFO] {fit_type} fit: All VAF peaks found in VAF peaks expected for ploidy {p}. Ploidy accepted')
            ploidy = p 
            break 
            
    else: # if cycle finished without finding ploidy 
        log_buffer.append(f'[WARNING] {fit_type} fit: VAF peaks do not match expected from ploides from 1, 2, 3. Contamination possible ')

    return True, ploidy, log_buffer 


def verdict_ploidy(first_data, first_ploidy, second_data, second_ploidy, pseudo_vector_in_sample):
    '''
    Combines results of fits to get final verdict 
    :param first_data: bool, if first fit had enough data for fit
    :param first_ploidy: int or bool, result of first fit 
    :param second_data: bool, if second fit had enough data for fit
    :param second_ploidy: int or bool, result of second fit 
    :param pseudo_vector_in_sample: float, persent of pseudogenic variants present in sample
    :return fit_ok: bool, if any fits had enough data
    :return ploidy_found: bool, if any fits were sucessfull 
    :return final_ploidy: int, ploidy if any
    :return verdict: str, text verdict 
    :return log_app_type: str, info or warning 
    '''
    
    pseudo_v_thres = 0.75
    fit_ok, ploidy_found, final_ploidy = (first_data | second_data), True, None
    
    if not fit_ok: # no data
        ploidy_found, log_app_type, final_ploidy = False, '[WARNING]', 2
        verdict = 'Not enough variants for sucessful fits. Contamination and ploidy will not be estimated. Default 2 is used'

    elif (not first_ploidy) & (not second_ploidy): # both ploidies not found 
        ploidy_found, log_app_type = False, '[WARNING]'
        verdict = 'All ploidy fits returned contamination or ploidy higher than 3.  Contamination risk'
    
    elif first_data & ~second_data: # only first fit present
        log_app_type, final_ploidy = '[INFO]', first_ploidy
        verdict = 'Ploidy fitted. No contamination risk'
        
    elif ~first_data & second_data: # only second fit present 
        log_app_type, final_ploidy = '[INFO]', second_ploidy
        verdict = 'Ploidy fitted. No contamination risk'
        
    elif (first_ploidy | second_ploidy) & ~(first_ploidy & second_ploidy): # one fit good, other no ploidy 
        ploidy_found, log_app_type, final_ploidy = False, '[WARNING]', max(first_ploidy, second_ploidy)
        verdict = 'One of the two fits returned contamination or ploidy higher than 3. Contamination risk'
        
    else: # 2 good fits 
        if first_ploidy == second_ploidy:  # results match 
            log_app_type, final_ploidy = '[INFO]', first_ploidy
            verdict = 'Ploidies in both fits match. No contamination'
            
        elif (first_ploidy != 3) & (second_ploidy != 3):
            if pseudo_vector_in_sample > pseudo_v_thres:
                log_app_type, final_ploidy = '[WARNING]', max(first_ploidy, second_ploidy)
                verdict = 'Ploidies in both fits do not match, both ploidies lower than 3. High pseudogenic content, low contamination risk is present'
            else: 
                log_app_type, final_ploidy = '[INFO]', max(first_ploidy, second_ploidy)
                verdict = 'Ploidies in both fits do not match, both ploidies lower than 3. No contamination risk'
        else: 
            if pseudo_vector_in_sample > pseudo_v_thres: 
                log_app_type, final_ploidy = '[WARNING]', max(first_ploidy, second_ploidy)
                verdict = 'Ploidies in both fits do not match, one of ploidies higher than 3. High pseudogenic content, contamination risk'
            else: 
                log_app_type, final_ploidy = '[INFO]', max(first_ploidy, second_ploidy)
                verdict = 'Ploidies in both fits do not match, one of ploidies higher than 3. No contamination risk. Ploidy detection may be unstable'

    return fit_ok, ploidy_found, final_ploidy, verdict, log_app_type


if __name__ == '__main__':  

    parser = argparse.ArgumentParser(description='Checks amplicon contamination, checks contamination with pseudogenic variants specifically')
    parser.add_argument('--hc_vcf', type=str, help='HaplotypeCaller vcf', required=True) # later may switch to mpileup results 
    parser.add_argument('--pseudovcf', type=str, help='Differences in amp sequence from pseudogene', required=False, default=None)
    parser.add_argument('--output_prefix', type=str, help='Output_prefix', required=True)
    args = parser.parse_args()
    print("PROGRAM ARGUMENTS: ", vars(args))

    vcf = read_vcf(args.hc_vcf)
    gt_vector = genotype_vectors(vcf)
    
    if args.pseudovcf is not None: 
        pseudo_vcf = read_vcf(link)
        pseudogene_vector = genotype_vectors(pseudo_vcf)
        intersection = pseudogene_vector.index.intersection(gt_vector.index)
        pseudo_vector_in_sample, nonpseudo_vector_in_sample = gt_vector.loc[intersection], gt_vector.drop(intersection)
        fit_type = 'Non-pseudogenic variants' 
    else: 
        pseudo_vector_in_sample, nonpseudo_vector_in_sample = [], gt_vector
        fit_type = 'All variants' 
    
    # do contamination check on all SNP or non-pseudogenic
    log_buffer = []
    first_data, first_ploidy, log_buffer = estimate_peak_fits(nonpseudo_vector_in_sample['VAF'], log_buffer, args.output_prefix, fit_type=fit_type)

    # # do contamination check on pseudo
    if len(pseudo_vector_in_sample) > 0: 
        second_data, second_ploidy, log_buffer = estimate_peak_fits(pseudo_vector_in_sample['VAF'], log_buffer, args.output_prefix, fit_type='Pseudogenic variants')
        percent_of_pseudogenic = len(pseudo_vector_in_sample)/len(pseudogene_vector)*100 
    else: 
        second_data, second_ploidy, percent_of_pseudogenic = False, False, None

    # if quality of fit was low because of split, try them all 
    if (len(pseudo_vector_in_sample) > 0) & ~first_data & ~second_data:
        first_data, first_ploidy, log_buffer = estimate_peak_fits(gt_vector['VAF'], log_buffer, args.output_prefix)

    verdict_results = {j: i for i, j in zip(verdict_ploidy(first_data, first_ploidy, second_data, second_ploidy, pseudo_vector_in_sample), 
                                            ['Enough_data_in_fit', 'Ploidy_well_fitted', 'ploidy', 'verdict', 'log_app_type'])}

    log_buffer.append(' '.join([verdict_results['log_app_type'], verdict_results['verdict']]))

    # write logs 
    with open(f'{args.output_prefix}.contamination_ploidy_check.log', 'w') as f: 
        for i in log_buffer:
            f.write(i + '\n')

    pd.Series(verdict_results).to_csv(f'{args.output_prefix}.contamination_ploidy_results.txt', sep='\t')