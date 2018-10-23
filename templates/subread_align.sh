#!/usr/bin/env bash

subread-align \
  -i ${idxmeta.target} \
  -r ${r1} \
  -R ${r2} \
  -t 0 \
  --SAMoutput \
  -T ${task.cpus} \
  > sam

  # subjunc \
  # -i ${idxmeta.target} \
  # -r ${r1} \
  # -R ${r2} \
  # -t 0 \
  # --SAMoutput \
  # --allJunctions \
  # -T ${task.cpus} \
  # > sam
