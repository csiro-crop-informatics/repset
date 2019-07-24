#!/usr/bin/env bash

bwa mem \
  -t ${task.cpus} \
  ${ALIGN_PARAMS} \
  ${idxmeta.target} \
  ${reads} \
  > out.sam
