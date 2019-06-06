#!/usr/bin/env bash

subread-align \
  -i ${idxmeta.target} \
  -r ${r1} \
  -R ${r2} \
  -t 0 \
  --SAMoutput \
  -T ${task.cpus} \
  ${ALIGN_PARAMS} \
  > sam

  # subjunc \
  # -i ${idxmeta.target} \
  # -r ${r1} \
  # -R ${r2} \
  # -t 0  (0 is RNA, 1 is DNA)
  # --SAMoutput \
  # --allJunctions \
  # -T ${task.cpus} \
  # > sam
