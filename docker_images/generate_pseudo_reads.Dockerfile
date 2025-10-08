FROM python:3.10

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

RUN curl -L https://github.com/ncsa/NEAT/archive/refs/tags/3.4.tar.gz | tar jxf -
RUN pip install -f NEAT-3.4/requirements.txt
