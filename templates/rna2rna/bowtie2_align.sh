#!/usr/bin/env bash

bowtie2 \
  -p ${task.cpus} \
  -x ${idxmeta.target} \
  -1 ${reads[0]} \
  -2 ${reads[1]} \
  --threads  ${task.cpus} \
  ${ALIGN_PARAMS} \
  > out.sam

