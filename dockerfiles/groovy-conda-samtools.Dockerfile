FROM groovy:3.0@sha256:ef1b4c220b71f85d264b742b7de2143a2cc824267756ad75aa818980da5797f5

ARG ANACONDA_VERSION=4.6.14

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends \
      libkeyutils-dev \
      ca-certificates \
      procps \
      libfontconfig1 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-${ANACONDA_VERSION}-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

ENV PATH /opt/conda/bin:$PATH

ARG SAMTOOLS_VERSION=1.9

RUN conda install --override-channels -c conda-forge -c bioconda -c default samtools=${SAMTOOLS_VERSION}