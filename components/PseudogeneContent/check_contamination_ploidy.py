import pandas as pd
import numpy as np 
import argparse 


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
    
    # add DP detection 
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
    
    return gt_vector


def extract_vaf_vector(gt_vector, gt, include_variants=[], exclude_variants=[]): 

    vaf_vector = {}
    if len(include_variants) > 0: 
        for variant_id in include_variants: 
            if (gt_vector[variant_id][0] == gt) & (variant_id not in exclude_variants): 
                vaf_vector[variant_id] = gt_vector[variant_id][1]
    else: 
        for variant_id in gt_vector: 
            if (gt_vector[variant_id][0] == gt) & (variant_id not in exclude_variants): 
                vaf_vector[variant_id] = gt_vector[variant_id][1]

    return vaf_vector 


def choose_best_gm(baf_input, max_ploidy=4): 
    '''
    Chooses best gaussian mixture for given baf input 
    :param baf_input: pd.Series, vector of bafs for sample
    :param max_ploidy: int, maximal ploidy allowed in fit 
    :return gm: sklearn.GaussianMixture, fitted model 
    '''
    baf_input = np.array(baf_input).reshape(-1, 1)
    
    # chose best gaussian
    bic, aic = {}, {}
    for i in range(1, max_ploidy+1): 
        gm = GaussianMixture(n_components=i, random_state=0).fit(baf_input) 
        bic[i] = gm.bic(baf_input)
        aic[i] = gm.aic(baf_input)
        
    # choose best components 
    best_n = (pd.Series(bic) + pd.Series(aic)).idxmin()
    
    print((pd.Series(bic) + pd.Series(aic)).sort_values())

    gm = GaussianMixture(n_components=best_n, random_state=0).fit(baf_input) 
    
    return best_n, gm 


def estimate_peak_fits(all_vaf_vector, homo_vaf_vector, hetero_vaf_vector): 
    best_n, gm = choose_best_gm(all_vaf_vector)
    


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
         # split pseudogenic check & do contamination
        pseudo_vcf = read_vcf(link)
        pseudogene_vector = genotype_vectors(pseudo_vcf)
    else: 
        # do contamination check on all SNP 
       