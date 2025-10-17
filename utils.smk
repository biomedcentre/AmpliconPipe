import glob
import os

# ===variables===

rep_location = config['repository']
reference_folder = os.path.dirname(config['references']['amplicon'])
reference_name = os.path.basename(config['references']['amplicon'])
reference_prefix = reference_name.split('.fasta')[0] 


pseudo_coords = config['references']['pseudogene_coordinates']
pseudo_vcf_use = f'--pseudovcf /pipeline/output/{reference_prefix}.vs.{pseudo_coords}.difference.vcf' if pseudo_coords is not None else '' 
multiqc_config = 'multiqc_config_template_pseudogenic.yaml' if pseudo_coords is not None else 'multiqc_config_template.yaml'
if config['references']['full_genome']:
    ref_name = os.path.basename(config['references']['full_genome'])
    full_ref_dir = os.path.dirname(config['references']['full_genome'])

bed_name = os.path.basename(config['references']['amplicon_bed'])
GID = config['user_settings']['GID']
output_folder = config['output']['output_folder']
threads = config['run_settings']['threads']
input_folder = config['input']['input_folder']


# ===functions===

def filename_to_sample(name):
    basename = os.path.basename(name)

    for i in config['input']['fastq_extensions']:
        if i in basename: 
            basename = basename.split(i)[0]
            break

    if '_L' in basename:
        # Find position of _L followed by 3 digits and _R1 or _R2
        import re
        pattern = r'_L\d{3}_R[12]'
        match = re.search(pattern, basename)
        if match:
            sample_name = basename[:match.start()]
            return sample_name

    if '_R1' in basename:
        sample_name = basename.split('_R1')[0]
        return sample_name
    elif '_R2' in basename:
        sample_name = basename.split('_R2')[0]
        return sample_name

    return basename 


def get_sample_names(input_folder): 
    '''
    Detects all sample names in input folder 
    :param: input_folder: str, given input folder
    '''
    all_files = []
    for i in config['input']['fastq_extensions']:
        all_files = all_files + glob.glob(os.path.join(input_folder, f'*{i}'))

    sample_names = [filename_to_sample(name) for name in all_files]

    # find sample names that have at least 2 fastq 
    sample_counts = {}
    for s in sample_names: 
        if s in sample_counts: 
            sample_counts[s] += 1
        else: 
            sample_counts[s] = 1

    final_samples = [s for s in sample_counts if sample_counts[s] > 1]

    return final_samples


def get_fastq_input(wildcards, input_folder, rgroup):
    sample = wildcards.sample
    all_files = []
    for i in config['input']['fastq_extensions']:
        all_files = all_files + glob.glob(os.path.join(input_folder, f'{sample}*{rgroup}*{i}'))

    return sorted(all_files)[0]


def get_alt_need(wildcards): 
     with checkpoints.PloidyContaminationQC.get(sample=wildcards.sample).output.qc_result.open() as file:
        line = file.readlines()[1].strip('\n').split('\t')
        if (line[-1] != "[CRITICAL]") & (line[-3] > 2):
            os.path.join(output_folder, f"{wildcards.sample}.alt.ploidy.phased.vcf")
        else: 
            return [] 