#!/usr/bin/env bash

minimap2 \
  -ax sr \
  -t ${task.cpus} \
  ${ALIGN_PARAMS} \
  ${idxmeta.target}.mmi \
  ${reads} \
  > out.sam
