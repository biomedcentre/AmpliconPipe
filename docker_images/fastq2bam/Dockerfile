FROM debian:bullseye-slim

ARG BWA_VERSION=2.2.1
ARG FASTP_VERSION=0.23.2

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        curl \
        samtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /pipeline/reference
RUN mkdir -p /pipeline/input
RUN mkdir -p /pipeline/output
RUN mkdir -p /pipeline/tools

WORKDIR /pipeline/tools

RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_VERSION}/bwa-mem2-${BWA_VERSION}_x64-linux.tar.bz2 | tar jxf - && \
    mv bwa-mem2-${BWA_VERSION}_x64-linux bwa-mem2
RUN curl --retry 5 --retry-delay 10 --retry-max-time 60 -L http://opengene.org/fastp/fastp.${FASTP_VERSION} > ./fastp && \
    chmod +x ./fastp