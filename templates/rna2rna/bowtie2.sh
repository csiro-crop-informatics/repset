#!/usr/bin/env bash

bowtie2 \
  -p ${task.cpus} \
  -x ${ref} \
  -1 ${reads[0]} \
  -2 ${reads[1]} \
  --threads  ${task.cpus} \
  ${ALIGN_PARAMS} \
  > out.sam

