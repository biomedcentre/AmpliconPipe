# AmpliconPipe
Pipeline for processing amplicons 
<img width="2309" height="634" alt="AmpliconPipe (1)" src="https://github.com/user-attachments/assets/d7b0bf8a-74ab-4495-913a-9bbfc779c886" />

## Running
Modify config.yaml for your launch. 
Change input_dir to your input dir with fastq and output_dir to your chosen output dir. Add sample patterns include_patterns if needed.

Modify reference and repository links with actual links for your user and system.

To get group id for your user (GID varaible in config)
```
id -g
```
To launch 
```
snakemake -p --cores 15 --default-resources tmpdir="/your/folder/tmp"
```
