import os 
import argparse 
import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
from math import sqrt
from sklearn.mixture import GaussianMixture
import subprocess 

# general constants 
minimal_vaf_vector_len = 4
pseudogenic_thres = 50
# constants for estimate peak fits 
minimal_peak_size = 3
allowed_std = 10
max_peaks = 4
max_allowed = 0.08

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
    sns.histplot(vaf_vector, bins=15, binrange=(0, 1), ax=ax)
    sns.kdeplot(vaf_vector, ax=ax)
    plot_name = fit_type.replace('_', ' ').title()
    ax.set_title(f'{plot_name} VAF')
    ax.set_xlim(-0.05, 1.05)
    plt.rcParams['figure.dpi'] = 150
    plt.savefig(f'{prefix}.{fit_type}.VAF_mqc.png')


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
    peak_dist_thres = 1/(max_ploidy + 2)
    
    good_peaks_stats.index = good_peaks_stats['name']  
    peak_groups, current_group, current_stats = [], [good_peaks_stats['name'].iloc[0]], good_peaks_stats.iloc[0]
    
    for idx, row in good_peaks_stats.iloc[1:].iterrows():
        # extract peak intervals for groups 
        interval = (current_stats['mean'] - 2*allowed_std*current_stats['std'], current_stats['mean'] + 2*allowed_std*current_stats['std'])
        interval_row = (row['mean'] - 2*allowed_std*row['std'], row['mean'] + 2*allowed_std*row['std'])
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
            new_good_peaks_stats.append([''.join([str(g) for g in group]),
                                         vaf_vector[current_group_samples].mean(),
                                         vaf_vector[current_group_samples].std(),
                                         current_group_samples.sum()])

    new_good_peaks_stats = pd.DataFrame(new_good_peaks_stats, columns=['name', 'mean', 'std', 'samples'])
    
    return new_good_peaks_stats
            

def estimate_peak_fits(vaf_vector, log_buffer, max_ploidy=3, fit_type='All variants'):
    '''
    Fits GM for VAF vector and checks list of ploidies to see if any fits. 
    :param vaf_vector: pd.Series, vector of vafs for variants
    :param log_buffer: list, list of logged string about ploidy fits
    :param max_ploidy: int, max ploidy that will be tested. Default is 3 for better peak resolution
    :param fit_type: str, description of current set of variants 
    :return: bool, if there was enough data for fit
    :return ploidy: bool or int, fit result
    :return log_buffer: list, updated log_buffer 
    '''
    # choose best GM and filter before merge 
    best_n, gm = choose_best_gm(vaf_vector, max_components=min(max_peaks, len(vaf_vector) - 1))
    pred_labels = pd.Series(gm.predict(np.array(vaf_vector).reshape(-1, 1)), index=vaf_vector.index)
    
    peaks_stats = []
    for i, samples in pred_labels.value_counts().items(): 
        mu, std = gm.means_[:,0][i], sqrt(gm.covariances_[:,0][:,0][i])
        peaks_stats.append([i, mu, std, samples]) 
    peaks_stats = pd.DataFrame(peaks_stats, columns=['name', 'mean', 'std', 'samples']).sort_values(by='mean', ascending=False)
    
    concated_peaks_stats = concat_close_baf_peaks(peaks_stats, pred_labels, vaf_vector, allowed_std, max_ploidy).sort_values(by='mean', ascending=False)
    good_peaks_stats = concated_peaks_stats[concated_peaks_stats['samples'] > minimal_peak_size]
    
    if len(good_peaks_stats) == 0: # no peaks with enough variants 
        log_buffer.append(f'[WARNING] {fit_type} fit: each peak in the best fit contains less than {minimal_peak_size} variants. Contamination and ploidy will not be estimated')
        return False, False, log_buffer 
 
    ploidy = False 
    for p in range(1, max_ploidy+1): # iterating through ploidies to see if any fits 
        expected_vafs = vafs_for_ploidy(p)
        
        if len(expected_vafs) < len(good_peaks_stats):
            log_buffer.append(f'[INFO] {fit_type} fit: Number of VAF peaks more than expected for ploidy {p}. Ploidy rejected')
            continue 
            
        peak_match = {}
        for idx, row in good_peaks_stats.iterrows(): # looking for matches for peaks of fit
            for match in expected_vafs[::-1]: # starting from larger vafs as stats table also starts from larger stats 
                allowed_distance = min(allowed_std*row['std'], max_allowed)
                if (match > (row['mean'] - allowed_distance)) & (match < (row['mean'] + allowed_distance)):  
                    peak_match[row['name']] = match 
                    expected_vafs = [i for i in expected_vafs if i != match]
                    break # match found, skip iteration over expected vafs
            else: 
                log_buffer.append(f'[INFO] {fit_type} fit: No match found for peak with mean {row["mean"]} and std {row["std"]} in VAF peaks expected for ploidy {p}. Ploidy rejected')
                break # further iteration not required, one mismatching peaks is enough
                
        if len(peak_match) == len(good_peaks_stats):
            log_buffer.append(f'[INFO] {fit_type} fit: All VAF peaks found in VAF peaks expected for ploidy {p}. Ploidy accepted')
            ploidy = p 
            break 
            
    else: # if cycle finished without finding ploidy 
        log_buffer.append(f'[WARNING] {fit_type} fit: VAF peaks do not match expected from ploides from 1, 2, 3. Contamination possible ')

    return True, ploidy, log_buffer 


def process_vector_fit(vector, log_buffer, fit_results, fit_type):
    '''
    Fits ploidy distribution, save results into dict structure or returns that data not sufficient for fit
    :param vector: pd.DataFrame, gt table for sample
    :param log_buffer: list, list of logged string about ploidy fits
    :param fit_results: dict of dicts with fit results
    :param fit_type: str, description of current set of variants 
    :return fit_results: dict of dicts with fit results
    :return log_buffer: list, list of logged string about ploidy fits
    '''

    fit_results[fit_type] = {} 
    if (len(vector) < minimal_vaf_vector_len):
        log_buffer.append(f'[WARNING] {fit_type} fit: Not enough variants (<{minimal_vaf_vector_len}) for VAF peak analysis. Contamination and ploidy will not be estimated')
        fit_results[fit_type]['data'], fit_results[fit_type]['ploidy'] = False, False
    else: 
        fit_results[fit_type]['data'], fit_results[fit_type]['ploidy'], log_buffer = estimate_peak_fits(vector['VAF'], log_buffer, fit_type=fit_type)

    return fit_results, log_buffer


def plot_vaf_distributions_for_fits(vectors, mode, prefix):
    '''
    Plots fits depending on the mode.
    :param vectors: dict of gt tables
    :param mode: str, mode
    :param prefix: str, path to output
    '''

    names_of_plots = {'all': 'all_variants', 'pseudo': 'pseudogenic_variants', 'nonpseudo': 'non_pseudogenic_variants'}
    if mode == 'all': 
        plot_vaf_vector(vectors['all']['VAF'], prefix, names_of_plots['all']) 
    else: 
        for k in vectors: 
            plot_vaf_vector(vectors[k]['VAF'], prefix, names_of_plots[k]) 


def ploidy_fits(vectors, mode, log_buffer):
    '''
    Fits ploidy using vaf vectors. Switches mode if detectes not sufficient data in fit, does additional fits if needeed 
    :param vectors: dict of gt tables
    :param mode: str, mode
    :param log_buffer: list, list of logged string about ploidy fits
    '''
    
    fit_results = {}
    if mode == 'all': # if all mode is chosen
        fit_results, log_buffer = process_vector_fit(vectors['all'], log_buffer, fit_results, 'all')
    else: 
        fit_results, log_buffer = process_vector_fit(vectors['pseudo'], log_buffer, fit_results, 'pseudo')
        fit_results, log_buffer = process_vector_fit(vectors['nonpseudo'], log_buffer, fit_results, 'nonpseudo')
        
        # check if additional all fit is needed 
        if ((not fit_results['pseudo']['data']) | (not fit_results['nonpseudo']['data'])): 
            mode = 'all'
            log_buffer.append('[INFO]: One of two fits returned insufficient data for getting confident peaks. Switching to all mode despite pseudogene sequence given') 
            fit_results, log_buffer = process_vector_fit(vectors['all'], log_buffer, fit_results, 'all')
        elif (fit_results['pseudo']['ploidy'] == 1) | (fit_results['nonpseudo']['ploidy'] == 1):
            log_buffer.append('[INFO]: At least one of two fits returned ploidy 1. Additionaly fitting all mode to ensure that no peaks got missed due to division of data')
            fit_results, log_buffer = process_vector_fit(vectors['all'], log_buffer, fit_results, 'all')

    return fit_results, log_buffer, mode


def pseudogenic_content_verdict(percent_of_pseudogenic, pseudogenic_thres=50, mode='all'):
    '''
    Asseses pseudogenic content if at least one of two fits returned contamination and forms verdict 
    :param percent_of_pseudogenic: float, percent of pseudogenic variants out of all possible pseudogenic variants
    :param pseudogenic_thres: float, percent of pseudogenic variants out of all possible pseudogenic variants
    :param mode: str, mode
    :return verdict: str, verdict line
    '''
    if mode != 'all':
        insert = 'One of two fits'
    else: 
        instert = 'Fit(s)'

    if percent_of_pseudogenic < pseudogenic_thres:
        verdict = f'{insert} returned contamination or ploidy higher than 3. Pseudogenic content low. Non-pseudogenic contamination risk'
    else:
        verdict = f'{insert} returned contamination or ploidy higher than 3. Pseudogenic content high. Pseudogenic contamination risk'

    return verdict 


def verdict_ploidy(fit_results, log_buffer, mode, percent_of_pseudogenic):
    '''
    :param fit_results: dict of dicts with fit results
    :param log_buffer: list, list of logged string about ploidy fits
    :param mode: str, mode
    :param percent_of_pseudogenic: float, percent of pseudogenic variants out of all possible pseudogenic variants
    '''
    fit_ok = True
    
    if mode == 'all':
        if not fit_results['all']['data']:
            fit_ok = False 
            ploidy_found, log_app_type, final_ploidy = False, '[WARNING]', 2
            verdict = 'Not enough variants for confident peaks in fit. Contamination and ploidy will not be estimated. Default 2 is used'
            
        elif not bool(fit_results['all']['ploidy']):
            ploidy_found, log_app_type, final_ploidy = False, '[CRITICAL]', None 
            verdict = pseudogenic_content_verdict(percent_of_pseudogenic)
        
        else:
            ploidy_found, log_app_type, final_ploidy = True, '[INFO]', fit_results['all']['ploidy']
            verdict = 'Ploidy fitted. No contamination risk'

    else: 
        if (not bool(fit_results['pseudo']['ploidy'])) & (not bool(fit_results['nonpseudo']['ploidy'])):
            ploidy_found, log_app_type, final_ploidy = False, '[CRITICAL]', None
            verdict = pseudogenic_content_verdict(percent_of_pseudogenic)
        
        elif (not bool(fit_results['pseudo']['ploidy'])) | (not bool(fit_results['nonpseudo']['ploidy'])):
            ploidy_found, log_app_type, final_ploidy = True, '[WARNING]', max(fit_results['pseudo']['ploidy'], fit_results['nonpseudo']['ploidy'])
            verdict = pseudogenic_content_verdict(percent_of_pseudogenic, mode=mode)

        elif fit_results['pseudo']['ploidy'] != fit_results['nonpseudo']['ploidy']:
            
            if percent_of_pseudogenic > pseudogenic_thres:
                ploidy_found, log_app_type, final_ploidy = True, '[WARNING]', max(fit_results['pseudo']['ploidy'], fit_results['nonpseudo']['ploidy'])
                verdict = 'Ploidies in both fits do not match, both ploidies equal or lower than 3. High pseudogenic content, low contamination risk is present'
            else: 
                ploidy_found, log_app_type, final_ploidy = True, '[INFO]', max(fit_results['pseudo']['ploidy'], fit_results['nonpseudo']['ploidy'])
                verdict = 'Ploidies in both fits do not match, both ploidies  equal or lower than 3. No contamination risk'

        elif (fit_results['pseudo']['ploidy'] == 1) | (fit_results['nonpseudo']['ploidy'] == 1):
            
            if fit_results['all']['data'] & bool(fit_results['all']['ploidy']):
                ploidy_found, log_app_type, final_ploidy = True, '[INFO]', fit_results['all']['ploidy']
                verdict = 'Ploidy 1 in both fits, using all variants fit to determine joint ploidy'
            else: 
                ploidy_found, log_app_type, final_ploidy = True, '[WARNING]', 1
                verdict = 'Ploidy 1 in both fits, but all variants fit failed. Ploidy may be unconfirmed'

        else: 
            ploidy_found, log_app_type, final_ploidy = True, '[INFO]', max(fit_results['pseudo']['ploidy'], fit_results['nonpseudo']['ploidy'])
            verdict = 'Ploidies in both match, both ploidies equal or lower than 3. No contamination risk'
        
    return fit_ok, ploidy_found, final_ploidy, verdict, log_app_type          


if __name__ == '__main__':  

    parser = argparse.ArgumentParser(description='Checks amplicon contamination, checks contamination with pseudogenic variants specifically')
    parser.add_argument('--hc_vcf', type=str, help='HaplotypeCaller vcf', required=True) # later may switch to mpileup results 
    parser.add_argument('--pseudovcf', type=str, help='Differences in amp sequence from pseudogene', required=False, default=None)
    parser.add_argument('--output_prefix', type=str, help='Output_prefix', required=True)
    parser.add_argument('--gid', type=str, default='0')
    parser.add_argument('--container', type=str, default='singularity')
    args = parser.parse_args()
    print("PROGRAM ARGUMENTS: ", vars(args))

    vcf = read_vcf(args.hc_vcf)
    vectors = {'all': genotype_vectors(vcf)}
    log_buffer = []

    mode, pseudo_percent = 'all', None
    if args.pseudovcf is not None: 
        pseudogene_vector = genotype_vectors(read_vcf(args.pseudovcf))
        intersection = pseudogene_vector.index.intersection(vectors['all'].index)
        vectors['pseudo'], vectors['nonpseudo'] = vectors['all'].loc[intersection], vectors['all'].drop(intersection)
        if (len(vectors['pseudo']) >= minimal_vaf_vector_len) & (len(vectors['nonpseudo']) >= minimal_vaf_vector_len):
            mode = 'separate'
        else: 
            log_buffer.append(f'[INFO]: One of two modes (pseudo/nonpseudo) has <{minimal_vaf_vector_len} variants. Switching to all mode despite pseudogene sequence given')
        percent_of_pseudogenic = len(vectors['pseudo'])/len(pseudogene_vector)*100 

    # fit ploidies and draw plots for fitted 
    fit_results, log_buffer, mode = ploidy_fits(vectors, mode, log_buffer)
    plot_vaf_distributions_for_fits(vectors, mode, args.output_prefix)

    verdict_results = {j: i for i, j in zip(verdict_ploidy(fit_results, log_buffer, mode, percent_of_pseudogenic), 
                                            ['Enough_data_in_fit', 'Ploidy_well_fitted', 'ploidy', 'verdict', 'log_app_type'])}
    verdict_results['mode'] = mode
    log_buffer.append(' '.join([verdict_results['log_app_type'], verdict_results['verdict']]))

    # write logs 
    with open(f'{args.output_prefix}.contamination_ploidy_check.log', 'w') as f: 
        for i in log_buffer:
            f.write(i + '\n')

    pd.DataFrame(verdict_results, index=[os.path.basename(args.output_prefix)]).to_csv(f'{args.output_prefix}.contamination_ploidy_results_mqc.txt', sep='\t')
    
    if args.container == 'docker':
        subprocess.run(['chown', '-Rc', f':{args.gid}', os.path.dirname(args.output_prefix)], capture_output=True, text=True, check=True)
    subprocess.run(['chmod', '-Rc', 'g+w,o-rwx', os.path.dirname(args.output_prefix)], capture_output=True, text=True, check=True)