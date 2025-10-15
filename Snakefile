import os 

configfile: "config.yaml"
include: "utils.smk"

samples = get_sample_names(input_folder)
print(samples)

rule all: 
    input: 
        expand(os.path.join(output_folder, "{sample}.contamination_ploidy_results.txt"), sample=samples),
        expand(os.path.join(output_folder, "{sample}.phased.vcf"), sample=samples)
        

rule fastq2bam:
    threads: int(config['run_settings']['threads'])
    input: 
        r1 = lambda wildcards: get_fastq_input(wildcards, input_folder, 'R1'),
        r2 = lambda wildcards: get_fastq_input(wildcards, input_folder, 'R2')
        ref_0123 = os.path.join(reference_folder, "{reference_name}.0123"),
        ref_amb = os.path.join(reference_folder, "{reference_name}.amb"),
        ref_ann = os.path.join(reference_folder, "{reference_name}.ann"),
        ref_bwt = os.path.join(reference_folder, "{reference_name}.bwt.2bit.64"),
        ref_pac = os.path.join(reference_folder, "{reference_name}.pac")
    output: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai"),
        fastp_json = os.path.join(output_folder, "{sample}.json")
    shell: 
        '''
        docker run --rm -v {input_folder}:{input_folder}:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} fastq2bam \
        components/fastq2bam/fastq2bam.sh /pipeline/reference/{reference_name} {input.r1} {input.r2} {threads}
        '''


rule HaplotypeCaller: 
    threads: int(config['run_settings']['threads'])
    input: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai"),
        ref_dict = os.path.join(reference_folder, "{reference_prefix}.dict")
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
        -v {reference_folder}:/pipeline/reference:ro \
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


rule index_reference: 
    threads: int(config['run_settings']['threads'])
    input:
        ref_fasta = config['references']['amplicon']
    output: 
        ref_0123 = os.path.join(reference_folder, "{reference_name}.0123"),
        ref_amb = os.path.join(reference_folder, "{reference_name}.amb"),
        ref_ann = os.path.join(reference_folder, "{reference_name}.ann"),
        ref_bwt = os.path.join(reference_folder, "{reference_name}.bwt.2bit.64"),
        ref_pac = os.path.join(reference_folder, "{reference_name}.pac")
    shell: 
        '''
        docker run --rm -v {reference_folder}:/pipeline/reference \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -e GID={GID} fastq2bam \
        components/index_bwa/index.sh /pipeline/reference/{reference_name} 
        '''


rule gatk_dict:
    threads: int(config['run_settings']['threads'])
    input:
        ref_fasta = config['references']['amplicon']
    output: 
        ref_dict = os.path.join(reference_folder, "{reference_prefix}.dict")
    shell:
        '''
        docker run --rm -v {reference_folder}:/pipeline/reference \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -e GID={GID} gatk \
        components/gatk_dict/gatk_dict.sh /pipeline/reference/{reference_name} /pipeline/reference/{reference_prefix}.dict
        '''


rule pseudogenic_synth_reads: 
    threads: int(config['run_settings']['threads'])
    input: 
        full_ref_dir = os.path.dirname(config['references']['full_genome']),
        ref_name = os.path.basename(config['references']['full_genome']),
    output: 
        pseudo_r1 = os.path.join(output_folder, "{pseudo_coords}.temp.pseudo_read1.fq.gz"),
        pseudo_r2 = os.path.join(output_folder, "{pseudo_coords}.temp.pseudo_read2.fq.gz")
    shell: 
        '''
        docker run --rm -v {input.full_ref_dir}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {output_folder}:/pipeline/output \
        -e GID=2009 -generate_pseudo_reads \
        components/GenePseudogeneDifference/generate_pseudo_reads_fastq.sh -P {input.pseudo_coords} \
        -R /pipeline/reference/{input.refname} -j {threads}
        '''


rule generate_pseudo_synth_bam: 
    threads: int(config['run_settings']['threads'])
    input: 
        pseudo_r1 = os.path.join(output_folder, "{pseudo_coords}.temp.pseudo_read1.fq.gz"),
        pseudo_r2 = os.path.join(output_folder, "{pseudo_coords}.temp.pseudo_read2.fq.gz"),
        ref_0123 = os.path.join(reference_folder, "{reference_name}.0123"),
        ref_amb = os.path.join(reference_folder, "{reference_name}.amb"),
        ref_ann = os.path.join(reference_folder, "{reference_name}.ann"),
        ref_bwt = os.path.join(reference_folder, "{reference_name}.bwt.2bit.64"),
        ref_pac = os.path.join(reference_folder, "{reference_name}.pac")
    output: 
        temp_bam = os.path.join(output_folder, "{pseudo_coords}.pseudoalign.bam") 
    shell:
        '''
        docker run --rm -v {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} fastq2bam \
        components/GenePseudogeneDifference/generate_pseudo_bam.sh -P {pseudo_coords} -j {threads} -a /pipeline/reference/{reference_name}
        '''


rule call_variants_different_in_gene: 
    threads: int(config['run_settings']['threads'])
    input:
        temp_bam = os.path.join(output_folder, "{pseudo_coords}.pseudoalign.bam") 
    output:
        
