FROM rsuchecki/miniconda3:4.6.14_050661b0ef92865fde5aea442f3440d1a7532659@sha256:c19f9684db9de4dd6852f942fa7a6a6e873d169c3bbe36ac3a365cbc458dee7e

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get install -y groovy

ARG SAMTOOLS_VERSION=1.9

RUN conda install --override-channels -c conda-forge -c bioconda -c default samtools=${SAMTOOLS_VERSION}
