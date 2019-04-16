FROM ubuntu:16.04
LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y wget unzip build-essential make libz-dev

#overwrite at build time with e.g. --build-arg VERSION=1.3.1
ARG VERSION=1.3.5

#Install dart
WORKDIR /dart
RUN wget https://github.com/rsuchecki/Dart/archive/v${VERSION}.tar.gz \
  && tar xzvf v${VERSION}.tar.gz \
  && cd Dart-${VERSION} \
  && make
ENV PATH "$PATH:/dart/Dart-${VERSION}"

# Cleanup
RUN rm /dart/v${VERSION}.tar.gz \
  && apt-get clean \
  && apt-get remove --yes --purge build-essential make libz-dev