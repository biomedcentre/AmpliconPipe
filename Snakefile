import os 

configfile: "config.yaml"
include: "utils.smk"

samples = get_sample_names(input_folder)
print(samples)

rule all: 
    input: 
        expand(os.path.join(output_folder, "{sample}.contamination_ploidy_results.txt"), sample=samples),
        expand(os.path.join(output_folder, "{sample}.phased.vcf"), sample=samples),
        expand(os.path.join(output_folder, "{sample}.multiqc.html"), sample=samples)
        

rule fastq2bam:
    threads: int(config['run_settings']['threads'])
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
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
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


checkpoint PloidyContaminationQC:
    input: 
        hc_vcf = os.path.join(output_folder, "{sample}.haplotype.caller.vcf.gz"),
        pseudo_diff_vcf = os.path.join(output_folder, f"{reference_prefix}.vs.{pseudo_coords}.difference.vcf") if pseudo_coords is not None else [] 
    output: 
        qc_result = os.path.join(output_folder, "{sample}.contamination_ploidy_results.txt")
    shell: 
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} python_and_whatshap \
        python3 components/PloidyContaminationQC/check_contamination_ploidy.py --hc_vcf /pipeline/input/{wildcards.sample}.haplotype.caller.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID} {pseudo_vcf_use}
        '''


rule index_reference: 
    threads: int(config['run_settings']['threads'])
    input:
        ancient(config['references']['amplicon'])
    output: 
        ref_0123 = os.path.join(reference_folder, f"{reference_name}.0123"),
        ref_amb = os.path.join(reference_folder, f"{reference_name}.amb"),
        ref_ann = os.path.join(reference_folder, f"{reference_name}.ann"),
        ref_bwt = os.path.join(reference_folder, f"{reference_name}.bwt.2bit.64"),
        ref_pac = os.path.join(reference_folder, f"{reference_name}.pac")
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
        ancient(config['references']['amplicon'])
    output: 
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
    shell:
        '''
        docker run --rm -v {reference_folder}:/pipeline/reference \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -e GID={GID} gatk \
        components/gatk_dict/gatk_dict.sh /pipeline/reference/{reference_name} /pipeline/reference/{reference_prefix}.dict
        '''


rule pseudogenic_synth_reads: 
    threads: int(config['run_settings']['threads'])
    output: 
        temp(os.path.join(output_folder, f"{pseudo_coords}.temp.pseudo_read1.fq.gz")),
        temp(os.path.join(output_folder, f"{pseudo_coords}.temp.pseudo_read2.fq.gz"))
    shell: 
        '''
        docker run --rm -v {full_ref_dir}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {output_folder}:/pipeline/output \
        -e GID=2009 generate_pseudo_reads \
        components/GenePseudogeneDifference/generate_pseudo_reads_fastq.sh -P {pseudo_coords} \
        -R /pipeline/input/{ref_name} -j {threads}
        '''


rule generate_pseudo_synth_bam: 
    threads: int(config['run_settings']['threads'])
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
        temp_bam = os.path.join(output_folder, f"{pseudo_coords}.pseudoalign.bam"), 
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
    output:
        pseudo_diff_vcf = os.path.join(output_folder, f"{reference_prefix}.vs.{pseudo_coords}.difference.vcf") 
    shell:
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} gatk \
        components/GenePseudogeneDifference/call_gatk_variants.sh -P {pseudo_coords} -j {threads} -a /pipeline/reference/{reference_name} 
        '''


rule HaplotypeCallerAltPloidy: 
    input: 
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai"),
        ref_dict = os.path.join(reference_folder, f"{reference_prefix}.dict")
        # ploidy = 3  might read it from file here
    output:
        hc_vcf_alt = os.path.join(output_folder, "{sample}.alt.ploidy.haplotype.caller.vcf.gz"),
        hc_vcf_alt_tbi = os.path.join(output_folder, "{sample}.alt.ploidy.haplotype.caller.vcf.gz.tbi")
    shell:
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} gatk \
        components/HaplotypeCallerAltPloidy/HC_alt.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads} 3
        '''


rule mpileupAltPloidy: 
    threads: int(config['run_settings']['threads'])
    input:
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
    output: 
        mpileup_alt = os.path.join(output_folder, "{sample}.alt.ploidy.bcftools.vcf.gz"),
        mpileup_alt_tbi = os.path.join(output_folder, "{sample}.alt.ploidy.bcftools.vcf.gz.tbi")
        # ploidy = 3  might read it from file here
    shell: 
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {reference_folder}:/pipeline/reference:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} mpileup \
        /pipeline/tools/components/mpileupAltPloidy/mpileup_altpl.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam {threads} 3
        '''


rule resolve_conflicts_AltPloidy: 
    input:
        hc_vcf_alt = os.path.join(output_folder, "{sample}.alt.ploidy.haplotype.caller.vcf.gz"),
        mpileup_alt = os.path.join(output_folder, "{sample}.alt.ploidy.bcftools.vcf.gz")
    output:
        final_vcf_alt = os.path.join(output_folder, "{sample}.alt.ploidy.conflict.resolved.vcf")
    shell: 
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {output_folder}:/pipeline/output \
        -e GID={GID} python_and_whatshap \
        python3 components/resolve_conflictsAltPloidy/resolve_conflitcs_alt.py --hc_vcf /pipeline/input/{wildcards.sample}.alt.ploidy.haplotype.caller.vcf.gz --pileup /pipeline/input/{wildcards.sample}.alt.ploidy.bcftools.vcf.gz --output_prefix /pipeline/output/{wildcards.sample} --gid {GID}
        '''


rule whatshap_polyphase: 
    threads: int(config['run_settings']['threads'])
    input: 
        final_vcf = os.path.join(output_folder, "{sample}.alt.ploidy.conflict.resolved.vcf"),
        bam = os.path.join(output_folder, "{sample}.sorted.dedup.bam"),
        bai = os.path.join(output_folder, "{sample}.sorted.dedup.bam.bai")
        # ploidy = 3  might read it from file here
    output:
        phased_vcf_alt = os.path.join(output_folder, "{sample}.alt.ploidy.phased.vcf")
    shell: 
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
        -v {rep_location}/components:/pipeline/tools/components:ro \
        -v {output_folder}:/pipeline/output \
        -v {reference_folder}:/pipeline/reference:ro \
        -e GID={GID} python_and_whatshap \
        bash components/whatshapAltPloidy/whatshap_polyphase.sh /pipeline/reference/{reference_name} /pipeline/input/{wildcards.sample}.sorted.dedup.bam  /pipeline/output/{wildcards.sample}.alt.ploidy.conflict.resolved.vcf {threads} 3
        '''


rule multiqc: 
    threads: int(config['run_settings']['threads'])
    input:
        fastp_json = os.path.join(output_folder, "{sample}.json"),
        regions = os.path.join(output_folder, "{sample}.regions.bed.gz"),
        thresholds = os.path.join(output_folder, "{sample}.thresholds.bed.gz"),
        qc_result = os.path.join(output_folder, "{sample}.contamination_ploidy_results.txt"),
        phased_vcf_alt = lambda wildcards: get_alt_need(wildcards) 
    output: 
        multi_qc_report = os.path.join(output_folder, "{sample}.multiqc.html")
    shell: 
        '''
        docker run --rm -v  {output_folder}:/pipeline/input:ro \
            -v {rep_location}/components:/pipeline/tools/components:ro \
            -v {output_folder}:/pipeline/output \
            -v {reference_folder}:/pipeline/reference:ro \
            -e GID={GID} python_and_whatshap \
            components/multiqc/multiqc.sh {wildcards.sample} /pipeline/tools/components/multiqc/{multiqc_config}
        '''