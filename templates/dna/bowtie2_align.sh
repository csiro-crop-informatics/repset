#!/usr/bin/env bash

bowtie2 \
  -p ${task.cpus} \
  -x ${idxmeta.target} \
  -1 1.fq.gz \
  -2 2.fq.gz \
  --threads  ${task.cpus} \
  ${ALIGN_PARAMS} \
  > out.sam

