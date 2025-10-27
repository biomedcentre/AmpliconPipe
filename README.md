# AmpliconPipe
Pipeline for processing amplicons 
<img width="4613" height="1517" alt="AmpliconPipe (3)" src="https://github.com/user-attachments/assets/8f0f2d39-522e-424f-aec4-8f34e26cb323" />

# Installation 

```
git clone 
cd AmpliconPipe
```
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

# Running
Modify config.yaml for your launch. 
Change input_dir to your input dir with fastq and output_dir to your chosen output dir. Add sample patterns include_patterns if needed.

Modify reference and repository links with actual links for your user and system.

To get group id for your user (GID varaible in config):
```
id -g
```
To launch 
```
snakemake -p --cores 45 --default-resources "tmpdir='/primary/data/zantysheva/projects/vdkn/AmpPipe_testing/tmp'" --configfile configs/config.yaml --restart-times 1
```
