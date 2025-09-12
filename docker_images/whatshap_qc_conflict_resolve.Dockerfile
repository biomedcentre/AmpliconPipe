FROM python:3.12

ENV VER_WHATSHAP="2.3"

RUN pip install numpy==1.25.0 pandas==2.0.3 scikit-learn==1.3.0 scipy==1.11.1 seaborn==0.13.2 matplotlib==3.8.3
RUN pip install --user whatshap==2.3 

RUN mkdir -p /pipeline/input
RUN mkdir -p /pipeline/output
RUN mkdir -p /pipeline/references
RUN mkdir -p /pipeline/tools