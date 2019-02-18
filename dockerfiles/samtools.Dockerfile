FROM rsuchecki/miniconda3:4.5.12_d42c6c234cbabb3737a145df6a52230cf2841923

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

ARG SAMTOOLS_VERSION=1.9

RUN conda install --override-channels -c conda-forge -c bioconda -c default samtools=${SAMTOOLS_VERSION}