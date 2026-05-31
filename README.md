# AmpliconPipe
Pipeline for processing long-range amplicons short-read NGS from fastq. Includes variant harmonization between different callers, copy number calling based on VAFs in amplicon and additional QC based on VAF. If copy number equal to 3 detected, calls and filters variants for this copy number. Supports either read-based phasing or family-based phasing logic, if ped file is supplied.
Supplied reference example is CYP21A2 long-range amplicon, but pipeline can be set to run on any reference via config.
<img width="1819" height="864" alt="AmpliconPipe" src="https://github.com/user-attachments/assets/7cb20dc6-81b6-42d2-8dd1-a3afdb734969" />

# Requirements 

Requirements are snakemake v7.32.4 or higher, make v4.4.1 and engine (docker v29.1.3 or singularity v3.8.7 or higher). Snakemake, make, and singularity can be installed with provided conda enviroment (see Installation). Docker cannot be installed with conda, see https://docs.docker.com/engine/install/ for installation. 

# Installation 

Clone repostory.
```
git clone https://github.com/skyincyan/AmpliconPipe.git
cd AmpliconPipe
```
If you want to use conda enviroment for snakemake and singularity: 
```
conda create env -f enviroment.yaml
```
Next, building images for engine.
If your chosen engine is singularity:
```
make singularity
```
Alternatively, if you will be using docker containers:
```
make docker 
```
Or to build both kinds of images: 
```
make
```

# Preparing Minimal Run

You need to modify config-template.yaml for your launch.

1) Change `input:input_folder` to your input dir with fastq and `output:output_folder` to your chosen output dir. Input and output dirs should not be the same. 

2) Add path to your reference `references:amplicon` and `references:amplicon_bed` reference bed (interval where variants will be calling).  

3) Set `user_settings:GID` group id variable so that files will be saved with your user as owner (very important if engine is docker). 
To get group id for your user run in command line:
```
id -g
```

4) Set `run_settings:engine` to docker or singularity, whichever you are going to use. 

# Extra Options in Preparing Run 

1) AmpliconPipe can utilize relatedness information for logical family phasing. Provide ped file path `input:ped_file`
   
2) You can choose if you want save only harmonized caller output, our all callers outputs for exploration later. `output:save_all_callers` can be set to True (save extra output) or False (do not save)
   
3) If your amplicon region had highly homologous region (like CYP21A2 has CYP21A1P or like SMN1 has SMN2), you can provide path to full genomic sequence `references:full_genome` and coordinates of this region in full genomic sequence `references:pseudogene_coordinates`. AmpliconPipe will utilize this information for extra presision in copy number calling and assesment if long-range amplicon was contaminated by product of said highly homologous regions. 

# Running

Launch example
```
snakemake -p --cores 15 --default-resources "tmpdir='/your/tmp/dir/tmp'" --configfile configs/your-run-config.yaml --restart-times 1
```

# Results 

__Files in output for a sample:__ 
File Postfix | If Always Present | Explanation 
------------- | ------------- | ------------- 
.sorted.dedup.bam | Always | Reads aligned to provided amplicon reference 
.mosdepth.global.dist.txt, .mosdepth.region.dist.txt, .mosdepth.summary.txt, .regions.bed.gz, .thresholds.bed.gz | Always | Mosdepth output (coverage stats) 
.json | Always | Fastp output 
.conflict.resolved.vcf.gz | Always | Harmonyzed variant calls; vcf combines HC, DV and bcftools records that were found to be true 
.variant_conflict_resolution.log | Always | Extended caller conflict resolution log 
.contamination_ploidy_results_mqc.txt | Always | Results of copy number calling and contamination check 
.contamination_ploidy_check.log | Always | Extended log of distribution fit for copy number calling and contamination check 
.VAF_mqc.png | Always | VAF distributions of pseudogenic/non-pseudogenic/all variants 
.phased.vcf.gz | Always | Read-based phased vcf (phased with whatshap) 
.multiqc.html | Always | Multiqc report, combining mosdepth, fastp, copy number calling, and contamination check 
.multiqc_data | Always | Folder with data for multiqc report 
.family.phased.vcf | If ped file was given, sample had family and is not copy number 3 | Vcf with family-based phasing results 
.raw.caller.output | If output of all callers is saved | Folder with output of all 3 callers 
.raw.caller.output/.bcftools.vcf.gz | If output of all callers is saved | Bcftools calls 
.raw.caller.output/.deepvariant.vcf | If output of all callers is saved | DeepVariant calls
.raw.caller.output/.visual_report.html | If output of all callers is saved | DeepVariant stats 
.raw.caller.output/.haplotype.caller.vcf.gz | If output of all callers is saved | HaplotypeCaller calls
.alt.ploidy.quality.checked.vcf.gz | If predicted copy number is 3 | Variant calls via HaplotypeCaller with ploidy 3 after additional filtering 
.alt.ploidy.quality.log  | If predicted copy number is 3 | Extended log of filtering HaplotypeCaller ploidy 3 calls 
.alt.ploidy.phased.vcf.gz | If predicted copy number is 3  | Read-based phased vcf (phased with whatshap) for ploidy 3 calls
.alt.ploidy.haplotype.caller.vcf.gz | If predicted copy number is 3 and if output of all callers is saved | Variant calls via HaplotypeCaller with ploidy 3
