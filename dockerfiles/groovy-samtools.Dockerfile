FROM groovy:3.0-jre-alpine@sha256:99356f5d9a525e64fccc7f84a992a970fbb210275c6689b7d75398ca01010d19

USER root

CMD ["/bin/sh"]

ENV SAMTOOLS_VERSION=1.9

RUN apk add --no-cache build-base zlib-dev bzip2-dev xz-dev ncurses-dev ca-certificates wget bash \
    && wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION}/ && make && cd .. \
    && mv samtools-${SAMTOOLS_VERSION}/samtools /bin/ \
    && rm -rf samtools* \ 
    && apk del build-base zlib-dev ca-certificates wget

USER groovy