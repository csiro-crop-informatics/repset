#!/usr/bin/env bash

bowtie2 \
  -p ${task.cpus} \
  -x ${idxmeta.target} \
  -1 ${r1} \
  -2 ${r2} \
  -f \
  --threads  ${task.cpus} \
  --local \
  > sam

