import os 

configfile: "config.yaml"
include: "utils.smk"

wildcard_constraints:
    sample = r'\w+'

samples = get_sample_names(input_folder)
print(samples)


rule all: 
    input: 
        expand(os.path.join(output_folder, "{sample}.contamination_ploidy_results_mqc.txt"), sample=samples),
        expand(os.path.join(output_folder, "{sample}.phased.vcf"), sample=samples),
        expand(os.path.join(output_folder, "{sample}.multiqc.html"), sample=samples)
        

rule fastq2bam:
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task =  int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input: 
        r1 = lambda wildcards: get_fastq_input(wildcards, input_folder, 'R1'),
        r2 = lambda wildcards: get_fastq_input(wildcards, input_folder, 'R2'),
        ref_0123 = os.path.join(reference_folder, f"{reference_name}.0123"),
        ref_amb = os.path.join(reference_folder, f"{reference_name}.amb"),
        ref_ann = os.path.join(reference_folder, f"{reference_name}.ann"),
        ref_bwt = os.path.join(reference_folder, f"{reference_name}.bwt.2bit.64"),
        ref_pac = os.path.join(reference_folder, f"{reference_name}.pac")
    output: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai"),
        fastp_json = os.path.join(output_folder, "{sample}.json")
    params: 
        container_name = get_full_container('fastq2bam')
    shell: 
        '''
        {container_run_command} \
        {bind} {input_folder}:{input_folder}:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        components/fastq2bam/fastq2bam.sh /pipeline/reference/{reference_name} {input.r1} {input.r2} {threads}
        '''


rule HaplotypeCaller: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    threads: int(config['run_settings']['threads'])
    input: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai"),
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
    output:
        hc_vcf = os.path.join(output_folder, "{sample}.haplotype.caller.vcf.gz"),
        hc_vcf_tbi = os.path.join(output_folder, "{sample}.haplotype.caller.vcf.gz.tbi")
    params: 
        container_name = get_full_container('gatk')
    shell: 
        '''
        {container_run_command} \
        {bind {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        components/HaplotypeCaller/HC.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads}
        '''


rule mosdepth:
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input:
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output: 
        regions = os.path.join(output_folder, "{sample}.regions.bed.gz"),
        thresholds = os.path.join(output_folder, "{sample}.thresholds.bed.gz")
    params: 
        container_name = get_full_container('mosdepth')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        /pipeline/tools/components/mosdepth/mosdepth.sh /pipeline/input/{wildcards.sample}.sorted.dedup.bam /pipeline/reference/{bed_name} {threads}
        '''


rule mpileup: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input:
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output: 
        mpileup = os.path.join(output_folder, "{sample}.bcftools.vcf.gz"),
        mpileup_tbi = os.path.join(output_folder, "{sample}.bcftools.vcf.gz.tbi")
    params: 
        container_name = get_full_container('mpileup')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        /pipeline/tools/components/mpileup/mpileup.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads}
        '''


rule DeepVariant:
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input:
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output: 
        dv_vcf = os.path.join(output_folder, "{sample}.deepvariant.vcf")
    params: 
        container_name = get_full_container('deepvariant')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {env} THREADS={threads} {params.container_name} \
        /pipeline/tools/components/DeepVariant/DV.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam /pipeline/reference/{bed_name}
        '''


rule resolve_conflicts: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input:
        hc_vcf = os.path.join(output_folder, "{sample}.haplotype.caller.vcf.gz"),
        mpileup = os.path.join(output_folder, "{sample}.bcftools.vcf.gz"),
        dv_vcf = os.path.join(output_folder, "{sample}.deepvariant.vcf")
    output:
        final_vcf = os.path.join(output_folder, "{sample}.conflict.resolved.vcf")
    params: 
        container_name = get_full_container('python_and_whatshap')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        python3 components/resolve_conflicts/resolve_conflicts.py --hc_vcf /pipeline/input/{wildcards.sample}.haplotype.caller.vcf.gz --dv_vcf /pipeline/input/{wildcards.sample}.deepvariant.vcf --pileup /pipeline/input/{wildcards.sample}.bcftools.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID}
        '''


rule whatshap: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input: 
        final_vcf = os.path.join(output_folder, "{sample}.conflict.resolved.vcf"),
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output:
        phased_vcf = os.path.join(output_folder, "{sample}.phased.vcf")
    params: 
        container_name = get_full_container('python_and_whatshap')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {output_folder}:/pipeline/output \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {env} GID={GID} {params.container_name} \
        bash components/whatshap/whatshap.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam  /pipeline/output/{wildcards.sample}.conflict.resolved.vcf
        '''


checkpoint PloidyContaminationQC:
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input: 
        hc_vcf = os.path.join(output_folder, "{sample}.haplotype.caller.vcf.gz"),
        pseudo_diff_vcf = os.path.join(output_folder, f"{reference_prefix}.vs.{pseudo_coords}.difference.vcf") if pseudo_coords is not None else [] 
    output: 
        qc_result = os.path.join(output_folder, "{sample}.contamination_ploidy_results_mqc.txt")
    params: 
        container_name = get_full_container('python_and_whatshap')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        python3 components/PloidyContaminationQC/check_contamination_ploidy.py --hc_vcf /pipeline/input/{wildcards.sample}.haplotype.caller.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID} {pseudo_vcf_use}
        '''


rule index_reference: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input:
        ancient(config['references']['amplicon'])
    output: 
        ref_0123 = os.path.join(reference_folder, f"{reference_name}.0123"),
        ref_amb = os.path.join(reference_folder, f"{reference_name}.amb"),
        ref_ann = os.path.join(reference_folder, f"{reference_name}.ann"),
        ref_bwt = os.path.join(reference_folder, f"{reference_name}.bwt.2bit.64"),
        ref_pac = os.path.join(reference_folder, f"{reference_name}.pac")
    params: 
        container_name = get_full_container('fastq2bam')
    shell: 
        '''
        {container_run_command} \
        {bind} {reference_folder}:/pipeline/reference \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {env} GID={GID} {params.container_name} \
        components/index_bwa/index.sh /pipeline/reference/{reference_name} 
        '''


rule gatk_dict:
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input:
        ancient(config['references']['amplicon'])
    output: 
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
    params: 
        container_name = get_full_container('gatk')
    shell: 
        '''
        {container_run_command} \
        {bind} {reference_folder}:/pipeline/reference \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {env} GID={GID} {params.container_name} \
        components/gatk_dict/gatk_dict.sh /pipeline/reference/{reference_name} /pipeline/reference/{reference_prefix}.dict
        '''


rule pseudogenic_synth_reads: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    output: 
        temp(os.path.join(output_folder, f"{pseudo_coords}.temp.pseudo_read1.fq.gz")),
        temp(os.path.join(output_folder, f"{pseudo_coords}.temp.pseudo_read2.fq.gz"))
    params: 
        container_name = get_full_container('generate_pseudo_reads')
    shell: 
        '''
        {container_run_command} \
        {bind} {full_ref_dir}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID=2009 {params.container_name} \
        components/GenePseudogeneDifference/generate_pseudo_reads_fastq.sh -P {pseudo_coords} \
        -R /pipeline/input/{ref_name} -j {threads}
        '''


rule generate_pseudo_synth_bam: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input: 
        pseudo_r1 = os.path.join(output_folder, f"{pseudo_coords}.temp.pseudo_read1.fq.gz"),
        pseudo_r2 = os.path.join(output_folder, f"{pseudo_coords}.temp.pseudo_read2.fq.gz"),
        ref_0123 = os.path.join(reference_folder, f"{reference_name}.0123"),
        ref_amb = os.path.join(reference_folder, f"{reference_name}.amb"),
        ref_ann = os.path.join(reference_folder, f"{reference_name}.ann"),
        ref_bwt = os.path.join(reference_folder, f"{reference_name}.bwt.2bit.64"),
        ref_pac = os.path.join(reference_folder, f"{reference_name}.pac")
    output: 
        temp(os.path.join(output_folder, f"{pseudo_coords}.pseudoalign.bam"))
    params: 
        container_name = get_full_container('fastq2bam')
    shell: 
        '''
        {container_run_command} \
        {bind} {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        components/GenePseudogeneDifference/generate_pseudo_bam.sh -P {pseudo_coords} -j {threads} -a /pipeline/reference/{reference_name}
        '''


rule call_variants_different_in_gene: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input:
        temp_bam = os.path.join(output_folder, f"{pseudo_coords}.pseudoalign.bam"), 
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
    output:
        pseudo_diff_vcf = os.path.join(output_folder, f"{reference_prefix}.vs.{pseudo_coords}.difference.vcf") 
   params: 
        container_name = get_full_container('gatk')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        components/GenePseudogeneDifference/call_gatk_variants.sh -P {pseudo_coords} -j {threads} -a /pipeline/reference/{reference_name} 
        '''


rule HaplotypeCallerAltPloidy: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai"),
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
    output:
        hc_vcf_alt = os.path.join(output_folder, "{sample}.alt.ploidy.haplotype.caller.vcf.gz"),
        hc_vcf_alt_tbi = os.path.join(output_folder, "{sample}.alt.ploidy.haplotype.caller.vcf.gz.tbi")
   params: 
        container_name = get_full_container('gatk')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        components/HaplotypeCallerAltPloidy/HC_alt.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads} 3
        '''


rule quality_checked_AltPloidy: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input:
        hc_vcf_alt = os.path.join(output_folder, "{sample}.alt.ploidy.haplotype.caller.vcf.gz"),
    output:
        final_vcf_alt = os.path.join(output_folder, "{sample}.alt.ploidy.quality.checked.vcf")
    params: 
        container_name = get_full_container('python_and_whatshap')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        python3 components/quality_checked_AltPloidy/quality_checked.py --hc_vcf /pipeline/input/{wildcards.sample}.alt.ploidy.haplotype.caller.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID}
        '''


rule whatshap_polyphase: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input: 
        final_vcf = os.path.join(output_folder, "{sample}.alt.ploidy.quality.checked.vcf"),
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
        # ploidy = 3  might read it from file here
    output:
        phased_vcf_alt = os.path.join(output_folder, "{sample}.alt.ploidy.phased.vcf")
    params: 
        container_name = get_full_container('python_and_whatshap')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {output_folder}:/pipeline/output \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {env} GID={GID} {params.container_name} \
        bash components/whatshapAltPloidy/whatshap_polyphase.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam  /pipeline/output/{wildcards.sample}.alt.ploidy.quality.checked.vcf {threads} 3
        '''


rule multiqc: 
    threads: int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings'][workflow.current_rule.name]['runtime']),
        mem_mb =  int(config['rules_settings'][workflow.current_rule.name]['mem_mb']),
        cpus_per_task = int(config['rules_settings'][workflow.current_rule.name]['cpus_per_task'])
    input:
        fastp_json = os.path.join(output_folder, "{sample}.json"),
        regions = os.path.join(output_folder, "{sample}.regions.bed.gz"),
        thresholds = os.path.join(output_folder, "{sample}.thresholds.bed.gz"),
        qc_result = os.path.join(output_folder, "{sample}.contamination_ploidy_results_mqc.txt"),
        phased_vcf_alt = lambda wildcards: get_alt_need(wildcards, output_folder) 
    output: 
        multi_qc_report = os.path.join(output_folder, "{sample}.multiqc.html")
   params: 
        container_name = get_full_container('python_and_whatshap')
    shell: 
        '''
        {container_run_command} \
            {bind} {rep_location}/components:/pipeline/tools/components:ro \
            {bind} {output_folder}:/pipeline/output \
            {env} GID={GID} {params.container_name} \
            components/multiqc/multiqc.sh {wildcards.sample} 
        '''