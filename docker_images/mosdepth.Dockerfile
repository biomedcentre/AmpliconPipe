FROM alpine:3.15

RUN apk add --no-cache \
    bash \
    wget

RUN mkdir -p /pipeline/input
RUN mkdir -p /pipeline/output
RUN mkdir -p /pipeline/tools
RUN mkdir -p /pipeline/regions

RUN wget https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VER}/mosdepth && \
    chmod +x /mosdepth

WORKDIR /pipeline/tools