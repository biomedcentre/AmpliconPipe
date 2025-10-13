ARG DV_VERSION=1.5.0

FROM google/deepvariant:${DV_VERSION}

RUN mkdir -p /pipeline/reference \
    mkdir -p /pipeline/input \
    mkdir -p /pipeline/output \
    mkdir -p /pipeline/tools

WORKDIR /pipeline/tools 

ENTRYPOINT ["bash"]