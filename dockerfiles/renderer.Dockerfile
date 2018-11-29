FROM rocker/r-ver:3.5.1
LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y pandoc pandoc-citeproc texlive texlive-latex-extra

RUN install2.r bookdown rticles rmarkdown