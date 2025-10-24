import glob
import re
import os

# ===variables===

# references
rep_location = os.path.dirname(workflow.snakefile)
reference_folder = os.path.dirname(config['references']['amplicon'])
reference_name = os.path.basename(config['references']['amplicon'])
reference_prefix = reference_name.split('.fasta')[0] 
bed_name = os.path.basename(config['references']['amplicon_bed'])


# options if pseudogene given 
pseudo_coords = config['references']['pseudogene_coordinates']
pseudo_vcf_use = f'--pseudovcf /pipeline/output/{reference_prefix}.vs.{pseudo_coords}.difference.vcf' if pseudo_coords is not None else '' 
if config['references']['full_genome']:
    ref_name = os.path.basename(config['references']['full_genome'])
    full_ref_dir = os.path.dirname(config['references']['full_genome'])


# input output 
input_folder = config['input']['input_folder']
output_folder = config['output']['output_folder']


# settings 
GID = config['user_settings']['GID']
save_all_callers = config['output']['save_all_callers']


# container settings 
container_type = config['run_settings']['engine']
container_run_command = 'docker run --rm' if config['run_settings']['engine'] == 'docker' else 'singularity run --writable-tmpfs'
bind = '-v' if config['run_settings']['engine'] == 'docker' else '-B'
env = '-e' if config['run_settings']['engine'] == 'docker' else '--env'


# ===functions===

def filename_to_sample(name):
    '''
    Extracts sample name from fastq filename 
    :param: name: str, fastq filename
    :return: sample_name: str, sample name 
    '''
    basename = os.path.basename(name)

    for i in config['input']['fastq_extensions']:
        if i in basename: 
            basename = basename.split(i)[0]
            break

    if '_L' in basename:
        # Find position of _L followed by 3 digits and _R1 or _R2
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


def find_pattern(name, patterns):
    '''
    Searches for patterns in the name, if any present, returns true 
    :param: name: str, name 
    :param: patterns: re patterns to fild 
    :return: bool
    '''
    for pattern in patterns: 
        match = re.search(pattern, name)
        if match:
            return True
    return False


def get_sample_names(input_folder): 
    '''
    Detects all sample names in input folder 
    :param: input_folder: str, given input folder
    :return: final_samples: list of str, samples
    '''
    all_files = []
    for i in config['input']['fastq_extensions']:
        all_files = all_files + glob.glob(os.path.join(input_folder, f'*{i}'))

    sample_names = [filename_to_sample(name) for name in all_files]

    # filter sample names for patterns 
    for p_type in ['include_patterns', 'exclude_patterns']:
        if len(config['input'][p_type]) > 0: 
            sample_names = [name for name in sample_names if find_pattern(name, config['input'][p_type])]

    # find sample names that have at least 2 fastq 
    sample_counts = {}
    for s in sample_names: 
        if s in sample_counts: 
            sample_counts[s] += 1
        else: 
            sample_counts[s] = 1

    unpaired = [s for s in sample_counts if sample_counts[s] == 1]
    if len(unpaired) > 0:
        # TO DO: add proper snakemake warning here 
        print('WARNING: ', unpaired, ' samples have only one fastq. Unpaired data will not be processed')

    final_samples = [s for s in sample_counts if sample_counts[s] > 1]

    return final_samples


def get_fastq_input(wildcards, input_folder, rgroup):
    '''
    Gets R1 or R2 fastq file for the sample 
    :param: wildcards: Snakemake wildcards
    :param: input_folder: str, folder of input data
    :param: rgroup: str, R1 or R2
    :return: str, fastq file
    '''
    sample = wildcards.sample
    all_files = []
    for i in config['input']['fastq_extensions']:
        all_files = all_files + glob.glob(os.path.join(input_folder, f'{sample}*{rgroup}*{i}'))

    return sorted(all_files)[0]


def get_alt_need(wildcards, output_folder): 
    '''
    Reads PloidyContaminationQC results and determines if alt ploidy launch is needed. If needed, returns path to add to rule input so that rules generated alt ploidy vcf can be launched. 
    :param: wildcards: Snakemake wildcards
    :param: output_folder: str, folder where output is written 
    :return: str or empty list, paths to add to rule input 
    '''
     with checkpoints.PloidyContaminationQC.get(sample=wildcards.sample).output.qc_result.open() as file:
        line = file.readlines()[1].strip('\n').split('\t')
        if (line[-1] != "[CRITICAL]"): 
            if int(line[-3]) > 2:
                return os.path.join(output_folder, f"{wildcards.sample}.alt.ploidy.phased.vcf.gz")
        
        return []


def get_full_container(container_name): 
    '''
    Gets full container name depending on engine type 
    :param: container_name: str, name of container
    :returns: str, full name 
    '''
    if config['run_settings']['engine'] == 'singularity': 
        return os.path.join(rep_location, 'singularity_images', f'{container_name}.sif')
    else:
        return f'amppipe/container_name'
