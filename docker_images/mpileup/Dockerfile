## Supply this value at build time: --build-arg threads=<integer>
ARG threads=4

FROM debian:bullseye-slim AS build

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


FROM debian:bullseye-slim AS production

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        tabix \
        libbz2-dev \
        libcurl4-gnutls-dev \
        liblzma-dev \
        zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /usr/local/bin/bcftools /usr/local/bin/bcftools

WORKDIR /

RUN mkdir -p /pipeline/input
RUN mkdir -p /pipeline/output
RUN mkdir -p /pipeline/tools
RUN mkdir -p /pipeline/reference