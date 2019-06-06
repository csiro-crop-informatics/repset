#!/usr/bin/env bash

biokanga align \
  --sfx ${idxmeta.target}.sfx \
  --in 1.fq.gz \
  --pair 2.fq.gz  \
  --out out.bam \
  --threads ${task.cpus} \
  ${ALIGN_PARAMS}