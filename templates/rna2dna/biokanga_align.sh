#!/usr/bin/env bash

biokanga align \
  --sfx ${idxmeta.target}.sfx \
  --in ${reads[0]} \
  --pair ${reads[1]}  \
  --out out.bam \
  --threads ${task.cpus} \
  ${ALIGN_PARAMS}