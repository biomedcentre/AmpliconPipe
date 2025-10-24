.PHONY: deepvariant_singularity fastq2bam_singularity gatk_singularity generate_pseudo_reads_singularity mosdepth_singularity mpileup_singularity python_and_whatshap_singularity deepvariant_docker fastq2bam_docker gatk_docker generate_pseudo_reads_docker mosdepth_docker mpileup_docker python_and_whatshap_docker

all: singularity docker
singularity: deepvariant_singularity fastq2bam_singularity gatk_singularity generate_pseudo_reads_singularity mosdepth_singularity mpileup_singularity python_and_whatshap_singularity
docker: deepvariant_docker fastq2bam_docker gatk_docker generate_pseudo_reads_docker mosdepth_docker mpileup_docker python_and_whatshap_docker


deepvariant_singularity:
	mkdir -p singularity_images
	singularity build --fakeroot singularity_images/deepvariant.sif singularity_recipies/deepvariant.def

fastq2bam_singularity:
	mkdir -p singularity_images
	singularity build --fakeroot singularity_images/fastq2bam.sif singularity_recipies/fastq2bam.def

gatk_singularity:
	mkdir -p singularity_images
	singularity build --fakeroot singularity_images/gatk.sif singularity_recipies/gatk.def

generate_pseudo_reads_singularity:
	mkdir -p singularity_images
	singularity build --fakeroot singularity_images/generate_pseudo_reads.sif singularity_recipies/generate_pseudo_reads.def

mosdepth_singularity:
	mkdir -p singularity_images
	singularity build --fakeroot singularity_images/mosdepth.sif singularity_recipies/mosdepth.def

mpileup_singularity:
	mkdir -p singularity_images
	singularity build --fakeroot singularity_images/mpileup.sif singularity_recipies/mpileup.def

python_and_whatshap_singularity:
	mkdir -p singularity_images
	singularity build --fakeroot singularity_images/python_and_whatshap.sif singularity_recipies/python_and_whatshap.def

deepvariant_docker:
	docker build docker_images/DeepVariant \
	-t amppipe/deepvariant

fastq2bam_docker:
	docker build docker_images/fastq2bam \
	-t amppipe/fastq2bam

gatk_docker:
	docker build docker_images/gatk \
	-t amppipe/gatk

generate_pseudo_reads_docker:
	docker build docker_images/generate_pseudo_reads \
	-t amppipe/generate_pseudo_reads

mosdepth_docker:
	docker build docker_images/mosdepth \
	-t amppipe/mosdepth

mpileup_docker:
	docker build docker_images/mpileup \
	-t amppipe/mpileup

python_and_whatshap_docker:
	docker build docker_images/python_and_whatshap \
	-t amppipe/python_and_whatshap