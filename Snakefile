import os

include: "utils.smk"

wildcard_constraints:
    sample = r'[\w-]+'

samples = get_sample_names(input_folder)
print(samples)
families = get_families(ped_file)

rule all: 
    input: 
        expand(os.path.join(output_folder, "{sample}.contamination_ploidy_results_mqc.txt"), sample=samples),
        expand(os.path.join(output_folder, "{sample}.phased.vcf.gz"), sample=samples),
        expand(os.path.join(output_folder, "{sample}.multiqc.html"), sample=samples), 
        expand(os.path.join(output_folder, "{family}.family.ped.txt"), family=families) if family_phasing else []
        

rule fastq2bam:
    threads: int(config['rules_settings']['fastq2bam']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['fastq2bam']['runtime']),
        mem_mb =  int(config['rules_settings']['fastq2bam']['mem_mb']),
        cpus_per_task =  int(config['rules_settings']['fastq2bam']['cpus_per_task'])
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
        container_name = get_full_container('fastq2bam'),
        directory_for_calling = os.path.join(output_folder, "{sample}.raw.caller.output")
    shell: 
        '''
        {container_run_command} \
        {bind} {input_folder}:{input_folder}:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        /pipeline/tools/components/fastq2bam/fastq2bam.sh /pipeline/reference/{reference_name} {input.r1} {input.r2} {threads} {container_type} {wildcards.sample}
        mkdir -p {params.directory_for_calling}
        '''


rule HaplotypeCaller: 
    threads: int(config['rules_settings']['HaplotypeCaller']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['HaplotypeCaller']['runtime']),
        mem_mb =  int(config['rules_settings']['HaplotypeCaller']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['HaplotypeCaller']['cpus_per_task'])
    input: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai"),
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
    output:
        hc_vcf = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.haplotype.caller.vcf.gz") if save_all_callers else temp(os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.haplotype.caller.vcf.gz")),
        hc_vcf_tbi = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.haplotype.caller.vcf.gz.tbi") if save_all_callers else temp(os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.haplotype.caller.vcf.gz.tbi"))
    params: 
        container_name = get_full_container('gatk')
    shell: 
        '''
        {container_run_command} \
        {bind} {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {reference_folder}:/pipeline/reference:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        /pipeline/tools/components/HaplotypeCaller/HC.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads} {container_type}
        '''


rule mosdepth:
    threads: int(config['rules_settings']['mosdepth']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['mosdepth']['runtime']),
        mem_mb =  int(config['rules_settings']['mosdepth']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['mosdepth']['cpus_per_task'])
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
        /pipeline/tools/components/mosdepth/mosdepth.sh /pipeline/input/{wildcards.sample}.sorted.dedup.bam /pipeline/reference/{bed_name} {threads} {container_type}
        '''


rule mpileup: 
    threads: int(config['rules_settings']['mpileup']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['mpileup']['runtime']),
        mem_mb =  int(config['rules_settings']['mpileup']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['mpileup']['cpus_per_task'])
    input:
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output: 
        mpileup = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.bcftools.vcf.gz") if save_all_callers else temp(os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.bcftools.vcf.gz")),
        mpileup_tbi = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.bcftools.vcf.gz.tbi") if save_all_callers else temp(os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.bcftools.vcf.gz.tbi"))
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
        /pipeline/tools/components/mpileup/mpileup.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads} {container_type}
        '''


rule DeepVariant:
    threads: int(config['rules_settings']['DeepVariant']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['DeepVariant']['runtime']),
        mem_mb =  int(config['rules_settings']['DeepVariant']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['DeepVariant']['cpus_per_task'])
    input:
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output: 
        dv_vcf = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.deepvariant.vcf") if save_all_callers else temp(os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.deepvariant.vcf"))
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
        /pipeline/tools/components/DeepVariant/DV.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam /pipeline/reference/{bed_name} {container_type}
        '''


rule resolve_conflicts: 
    threads: int(config['rules_settings']['resolve_conflicts']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['resolve_conflicts']['runtime']),
        mem_mb =  int(config['rules_settings']['resolve_conflicts']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['resolve_conflicts']['cpus_per_task'])
    input:
        hc_vcf = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.haplotype.caller.vcf.gz"),
        mpileup = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.bcftools.vcf.gz"),
        dv_vcf = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.deepvariant.vcf")
    output:
        final_vcf = temp(os.path.join(output_folder, "{sample}.conflict.resolved.vcf"))
    params: 
        container_name = get_full_container('python_and_whatshap'),
        input_prefix_in_container = "/pipeline/input/{sample}.raw.caller.output/{sample}"
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        python3 /pipeline/tools/components/resolve_conflicts/resolve_conflicts.py --hc_vcf {params.input_prefix_in_container}.haplotype.caller.vcf.gz --dv_vcf {params.input_prefix_in_container}.deepvariant.vcf --pileup {params.input_prefix_in_container}.bcftools.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID} --container {container_type}
        '''


rule whatshap: 
    threads: int(config['rules_settings']['whatshap']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['whatshap']['runtime']),
        mem_mb =  int(config['rules_settings']['whatshap']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['whatshap']['cpus_per_task'])
    input: 
        final_vcf = os.path.join(output_folder, "{sample}.conflict.resolved.vcf"),
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output:
        phased_vcf = os.path.join(output_folder, "{sample}.phased.vcf.gz"),
        phased_vcf_tbi = os.path.join(output_folder, "{sample}.phased.vcf.gz.tbi")
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
        bash /pipeline/tools/components/whatshap/whatshap.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam  /pipeline/output/{wildcards.sample}.conflict.resolved.vcf {container_type}
        '''


checkpoint PloidyContaminationQC:
    threads: int(config['rules_settings']['PloidyContaminationQC']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['PloidyContaminationQC']['runtime']),
        mem_mb =  int(config['rules_settings']['PloidyContaminationQC']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['PloidyContaminationQC']['cpus_per_task'])
    input: 
        hc_vcf = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.haplotype.caller.vcf.gz"),
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
        python3 /pipeline/tools/components/PloidyContaminationQC/check_contamination_ploidy.py --hc_vcf /pipeline/input/{wildcards.sample}.raw.caller.output/{wildcards.sample}.haplotype.caller.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID} {pseudo_vcf_use} --container {container_type}
        '''


rule index_reference: 
    threads: int(config['rules_settings']['index_reference']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['index_reference']['runtime']),
        mem_mb =  int(config['rules_settings']['index_reference']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['index_reference']['cpus_per_task'])
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
        /pipeline/tools/components/index_bwa/index.sh /pipeline/reference/{reference_name} {container_type}
        '''


rule gatk_dict:
    threads: int(config['rules_settings']['gatk_dict']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['gatk_dict']['runtime']),
        mem_mb =  int(config['rules_settings']['gatk_dict']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['gatk_dict']['cpus_per_task'])
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
        /pipeline/tools/components/gatk_dict/gatk_dict.sh /pipeline/reference/{reference_name} /pipeline/reference/{reference_prefix}.dict {container_type}
        '''


rule pseudogenic_synth_reads: 
    threads: int(config['rules_settings']['pseudogenic_synth_reads']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['pseudogenic_synth_reads']['runtime']),
        mem_mb =  int(config['rules_settings']['pseudogenic_synth_reads']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['pseudogenic_synth_reads']['cpus_per_task'])
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
        /pipeline/tools/components/GenePseudogeneDifference/generate_pseudo_reads_fastq.sh -P {pseudo_coords} \
        -R /pipeline/input/{ref_name} -j {threads} -c {container_type}
        '''


rule generate_pseudo_synth_bam: 
    threads: int(config['rules_settings']['generate_pseudo_synth_bam']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['generate_pseudo_synth_bam']['runtime']),
        mem_mb =  int(config['rules_settings']['generate_pseudo_synth_bam']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['generate_pseudo_synth_bam']['cpus_per_task'])
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
        /pipeline/tools/components/GenePseudogeneDifference/generate_pseudo_bam.sh -P {pseudo_coords} -j {threads} -a /pipeline/reference/{reference_name} -c {container_type}
        '''


rule call_variants_different_in_gene: 
    threads: int(config['rules_settings']['call_variants_different_in_gene']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['call_variants_different_in_gene']['runtime']),
        mem_mb =  int(config['rules_settings']['call_variants_different_in_gene']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['call_variants_different_in_gene']['cpus_per_task'])
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
        /pipeline/tools/components/GenePseudogeneDifference/call_gatk_variants.sh -P {pseudo_coords} -j {threads} -a /pipeline/reference/{reference_name} -c {container_type}
        '''


rule HaplotypeCallerAltPloidy: 
    threads: int(config['rules_settings']['HaplotypeCallerAltPloidy']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['HaplotypeCallerAltPloidy']['runtime']),
        mem_mb =  int(config['rules_settings']['HaplotypeCallerAltPloidy']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['HaplotypeCallerAltPloidy']['cpus_per_task'])
    input: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai"),
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
    output:
        hc_vcf_alt = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.alt.ploidy.haplotype.caller.vcf.gz") if save_all_callers else temp(os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.alt.ploidy.haplotype.caller.vcf.gz")),
        hc_vcf_alt_tbi = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.alt.ploidy.haplotype.caller.vcf.gz.tbi") if save_all_callers else temp(os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.alt.ploidy.haplotype.caller.vcf.gz.tbi"))
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
        /pipeline/tools/components/HaplotypeCallerAltPloidy/HC_alt.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads} 3 {container_type}
        '''


rule quality_checked_AltPloidy: 
    threads: int(config['rules_settings']['quality_checked_AltPloidy']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['quality_checked_AltPloidy']['runtime']),
        mem_mb =  int(config['rules_settings']['quality_checked_AltPloidy']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['quality_checked_AltPloidy']['cpus_per_task'])
    input:
        hc_vcf_alt = os.path.join(output_folder, "{sample}.raw.caller.output", "{sample}.alt.ploidy.haplotype.caller.vcf.gz"),
    output:
        final_vcf_alt = temp(os.path.join(output_folder, "{sample}.alt.ploidy.quality.checked.vcf"))
    params: 
        container_name = get_full_container('python_and_whatshap')
    shell: 
        '''
        {container_run_command} \
        {bind}  {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {output_folder}:/pipeline/output \
        {env} GID={GID} {params.container_name} \
        python3 /pipeline/tools/components/quality_checked_AltPloidy/quality_checked.py --hc_vcf /pipeline/input/{wildcards.sample}.raw.caller.output/{wildcards.sample}.alt.ploidy.haplotype.caller.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID} --container {container_type}
        '''


rule whatshap_polyphase: 
    threads: int(config['rules_settings']['whatshap_polyphase']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['whatshap_polyphase']['runtime']),
        mem_mb =  int(config['rules_settings']['whatshap_polyphase']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['whatshap_polyphase']['cpus_per_task'])
    input: 
        final_vcf = os.path.join(output_folder, "{sample}.alt.ploidy.quality.checked.vcf"),
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
        # ploidy = 3  might read it from file here
    output:
        phased_vcf_alt = os.path.join(output_folder, "{sample}.alt.ploidy.phased.vcf.gz"),
        phased_vcf_alt_tbi = os.path.join(output_folder, "{sample}.alt.ploidy.phased.vcf.gz.tbi")
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
        bash /pipeline/tools/components/whatshapAltPloidy/whatshap_polyphase.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam  /pipeline/output/{wildcards.sample}.alt.ploidy.quality.checked.vcf {threads} 3 --container {container_type}
        '''


rule relatedness_phasing: 
    threads: int(config['rules_settings']['relatedness_phasing']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['relatedness_phasing']['runtime']),
        mem_mb =  int(config['rules_settings']['relatedness_phasing']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['relatedness_phasing']['cpus_per_task'])
    input:
        expand(os.path.join(output_folder, "{sample_f}.contamination_ploidy_results_mqc.txt"), sample_f=samples_in_family(wildcards, ped_file)),
        expand(os.path.join(output_folder, "{sample_f}.conflict.resolved.vcf"), sample_f=samples_in_family(wildcards, ped_file)),
    output:
        temp(os.path.join(output_folder, "{family}.family.ped.txt")),
        expand(os.path.join(output_folder, "{sample_f}.family.phased.vcf"), sample_f=samples_in_family(wildcards, ped_file))
    params: 
        container_name = get_full_container('python_and_whatshap')
    shell: 
        '''
        {container_run_command} \
        {bind} {output_folder}:/pipeline/input:ro \
        {bind} {rep_location}/components:/pipeline/tools/components:ro \
        {bind} {output_folder}:/pipeline/output \
        {bind} {ped_folder}:/pipeline/reference:ro \
        {env} GID={GID} {params.container_name} \
        bash /pipeline/tools/components/phasing_relatedness/phasing_related.py --input_folder /pipeline/input --family {wildcards.family} \
        --ped_file /pipeline/reference/{ped_name} --output_prefix /pipeline/output --container {container_type}
        '''
        

rule multiqc: 
    threads: int(config['rules_settings']['multiqc']['cpus_per_task'])
    resources: 
        runtime =  int(config['rules_settings']['multiqc']['runtime']),
        mem_mb =  int(config['rules_settings']['multiqc']['mem_mb']),
        cpus_per_task = int(config['rules_settings']['multiqc']['cpus_per_task'])
    input:
        fastp_json = os.path.join(output_folder, "{sample}.json"),
        regions = os.path.join(output_folder, "{sample}.regions.bed.gz"),
        thresholds = os.path.join(output_folder, "{sample}.thresholds.bed.gz"),
        qc_result = os.path.join(output_folder, "{sample}.contamination_ploidy_results_mqc.txt"),
        phased_vcf = os.path.join(output_folder, "{sample}.phased.vcf.gz"),
        phased_vcf_alt = lambda wildcards: get_alt_need(wildcards, output_folder)
    output: 
        multi_qc_report = os.path.join(output_folder, "{sample}.multiqc.html")
    params: 
        container_name = get_full_container('python_and_whatshap'),
        directory_for_calling = os.path.join(output_folder, "{sample}.raw.caller.output")
    shell: 
        '''
        {container_run_command} \
            {bind} {rep_location}/components:/pipeline/tools/components:ro \
            {bind} {output_folder}:/pipeline/output \
            {env} GID={GID} {params.container_name} \
            /pipeline/tools/components/multiqc/multiqc.sh {wildcards.sample} {container_type}

        if [ {save_all_callers} = 'False' ]; then
            rm -rd {params.directory_for_calling}
        fi 
        '''