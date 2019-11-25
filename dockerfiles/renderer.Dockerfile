FROM rocker/tidyverse:3.6.1

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update \
  && apt-get -qq -y install --no-install-recommends pandoc pandoc-citeproc texlive texlive-latex-extra procps wget xzdec libssl-dev libxml2-dev \
  && apt-get -qq -y autoremove \
  && apt-get autoclean \
  && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN install2.r bookdown rticles rmarkdown kableExtra