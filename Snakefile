import os 

configfile: "config.yaml"

rep_location = config['repository']
reference_folder = os.path.dirname(config['references']['amplicon'])
reference_name = os.path.basename(config['references']['amplicon'])
bed_name = os.path.basename(config['references']['amplicon_bed'])
GID = config['user_settings']['GID']
output_folder = config['output']['output_folder']
threads = config['run_settings']['threads']
input_folder = config['input']['input_folder']

samples = ["770924470801_CYP21A2_S40"]

rule all: 
    input: 
        expand(os.path.join(output_folder, "{sample}.contamination_ploidy_results.txt"), sample=samples),
        expand(os.path.join(output_folder, "{sample}.phased.vcf"), sample=samples)
        

rule fastq2bam:
    threads: int(config['run_settings']['threads'])
    input: 
        r1 = os.path.join(input_folder, "{sample}_R1_001.fastq.gz"),
        r2 = os.path.join(input_folder, "{sample}_R2_001.fastq.gz")
    output: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai"),
        fastp_json = os.path.join(output_folder, "{sample}.json")
    shell: 
        '''
        docker run --rm -v /primary/data/zantysheva/projects/vdkn/AmpPipe_testing/in:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} fastq2bam \
        components/fastq2bam/fastq2bam.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}_R1_001.fastq.gz /pipeline/input/{wildcards.sample}_R2_001.fastq.gz {threads}
        '''


rule HaplotypeCaller: 
    threads: int(config['run_settings']['threads'])
    input: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output:
        hc_vcf = os.path.join(output_folder, "{sample}.haplotype.caller.vcf.gz"),
        hc_vcf_tbi = os.path.join(output_folder, "{sample}.haplotype.caller.vcf.gz.tbi")
    shell:
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} gatk \
        components/HaplotypeCaller/HC.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads}
        '''


rule mosdepth:
    threads: int(config['run_settings']['threads'])
    input:
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output: 
        regions = os.path.join(output_folder, "{sample}.regions.bed.gz"),
        thresholds = os.path.join(output_folder, "{sample}.thresholds.bed.gz")
    shell:
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} gatk \
        /pipeline/tools/components/mosdepth/mosdepth.sh /pipeline/input/{wildcards.sample}.sorted.dedup.bam /pipeline/reference/{bed_name} {threads}
        '''


rule mpileup: 
    threads: int(config['run_settings']['threads'])
    input:
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output: 
        mpileup = os.path.join(output_folder, "{sample}.bcftools.vcf.gz"),
        mpileup_tbi = os.path.join(output_folder, "{sample}.bcftools.vcf.gz.tbi")
    shell: 
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} mpileup \
        /pipeline/tools/components/mpileup/mpileup.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads}
        '''


rule DeepVariant:
    threads: int(config['run_settings']['threads'])
    input:
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output: 
        dv_vcf = os.path.join(output_folder, "{sample}.deepvariant.vcf")
    shell:
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} -e THREADS={threads} deepvariant \
        /pipeline/tools/components/DeepVariant/DV.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam /pipeline/reference/{bed_name}
        '''


rule resolve_conflicts: 
    input:
        hc_vcf = os.path.join(output_folder, "{sample}.haplotype.caller.vcf.gz"),
        mpileup = os.path.join(output_folder, "{sample}.bcftools.vcf.gz"),
        dv_vcf = os.path.join(output_folder, "{sample}.deepvariant.vcf")
    output:
        final_vcf = os.path.join(output_folder, "{sample}.conflict.resolved.vcf")
    shell: 
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} python_and_whatshap \
        python3 components/resolve_conflicts/resolve_conflicts.py --hc_vcf /pipeline/input/{wildcards.sample}.haplotype.caller.vcf.gz --dv_vcf /pipeline/input/{wildcards.sample}.deepvariant.vcf --pileup /pipeline/input/{wildcards.sample}.bcftools.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID}
        '''


rule whatshap: 
    input: 
        final_vcf = os.path.join(output_folder, "{sample}.conflict.resolved.vcf"),
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output:
        phased_vcf = os.path.join(output_folder, "{sample}.phased.vcf")
    shell: 
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} python_and_whatshap \
        bash components/whatshap/whatshap.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam  /pipeline/output/{wildcards.sample}.conflict.resolved.vcf
        '''


rule PloidyContaminationQC:
    input: 
        hc_vcf = os.path.join(output_folder, "{sample}.haplotype.caller.vcf.gz")
    output: 
        qc_result = os.path.join(output_folder, "{sample}.contamination_ploidy_results.txt")
    shell: 
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} python_and_whatshap \
        python3 components/PloidyContaminationQC/check_contamination_ploidy.py --hc_vcf /pipeline/input/{wildcards.sample}.haplotype.caller.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID}
        '''