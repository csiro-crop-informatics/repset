FROM rocker/r-ver:3.5.2

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends procps \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN R -e "install.packages('dplyr', repos = 'https://cran.csiro.au')"
RUN R -e "install.packages('readr', repos = 'https://cran.csiro.au')"
RUN R -e "install.packages('ggplot2', repos = 'https://cran.csiro.au')"
RUN R -e "install.packages('viridis', repos = 'https://cran.csiro.au')"
RUN R -e "install.packages('ggrepel', repos = 'https://cran.csiro.au')"
RUN R -e "install.packages('purrr', repos = 'https://cran.csiro.au')"
RUN R -e "install.packages('tidyr', repos = 'https://cran.csiro.au')"

#RUN install2.r pkg1 pgk2 pkg3 ...

