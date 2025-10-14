FROM python:3.12-trixie

RUN apt-get update && apt-get install -y \
    ca-certificates \   
    build-essential \
    tabix && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /pipeline/input
RUN mkdir -p /pipeline/output
RUN mkdir -p /pipeline/references
RUN mkdir -p /pipeline/tools

COPY requirements.txt /pipeline/tools

RUN pip install --upgrade pip setuptools wheel 
RUN pip install -r /pipeline/tools/requirements.txt

WORKDIR /pipeline/tools

ENTRYPOINT [""]