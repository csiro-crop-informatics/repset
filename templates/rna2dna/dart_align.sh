#!/usr/bin/env bash

dart -i ${idxmeta.target} \
  -f 1.fq.gz \
  -f2 2.fq.gz \
  -t ${task.cpus} \
  ${ALIGN_PARAMS} \
  -o out.sam