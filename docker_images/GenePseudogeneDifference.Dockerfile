FROM broadinstitute/gatk:4.5.0.0 AS build

## Supply this value at build time: --build-arg threads=<integer>
ARG threads

## https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#sort-multi-line-arguments
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        autoconf \
        gcc \
        git \
        libbz2-dev \
        libcurl4-gnutls-dev  \
        liblzma-dev \
        make \
        zlib1g-dev && \
    apt-get install -y --reinstall ca-certificates

RUN git clone --recurse-submodules https://github.com/samtools/htslib.git
RUN git clone https://github.com/samtools/bcftools.git

WORKDIR /bcftools

RUN autoheader && autoconf && ./configure
RUN make -j $threads && make install



FROM broadinstitute/gatk:4.5.0.0

RUN apt-get update --allow-insecure-repositories && \
    apt-get install -y --no-install-recommends \
        parallel && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /pipeline/reference
RUN mkdir -p /pipeline/input
RUN mkdir -p /pipeline/output
RUN mkdir -p /pipeline/tools

WORKDIR /pipeline/tools

RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_VERSION}/bwa-mem2-${BWA_VERSION}_x64-linux.tar.bz2 | tar jxf -
RUN curl -L https://github.com/ncsa/NEAT/archive/refs/tags/3.4.tar.gz | tar jxf -

