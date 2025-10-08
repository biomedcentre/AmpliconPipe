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