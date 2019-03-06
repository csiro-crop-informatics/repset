FROM rsuchecki/miniconda3:4.5.12_d42c6c234cbabb3737a145df6a52230cf2841923

RUN conda install --override-channels -c conda-forge -c bioconda -c default rnftools