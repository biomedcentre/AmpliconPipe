# AmpliconPipe
Pipeline for processing long-range amplicons short-read NGS from fastq. Includes variant harmonization between different callers, copy number calling based on VAFs in amplicon and additional QC based on VAF. If copy number equal to 3 detected, calls and filters variants for this copy number. Supports either read-based phasing or family-based phasing logic, if ped file is supplied.
Supplied reference example is CYP21A2 long-range amplicon, but pipeline can be set to run on any reference via config.
<img width="4613" height="1517" alt="AmpliconPipe (3)" src="https://github.com/user-attachments/assets/8f0f2d39-522e-424f-aec4-8f34e26cb323" />

# Requirements 

Requirements are snakemake v7.32.4 and engine (docker v29.1.3 or singularity v3.8.7). Snakemake and singularity can be installed with provided conda enviroment (see Installation). Docker cannot be installed with conda, see https://docs.docker.com/engine/install/ for installation. 

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