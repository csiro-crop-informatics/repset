FROM rsuchecki/miniconda3:4.5.12_35a77a74877e5577f3c1110e3a6773d85c5c10c1

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

ARG SAMTOOLS_VERSION=1.9

RUN conda install --override-channels -c conda-forge -c bioconda -c default samtools=${SAMTOOLS_VERSION}