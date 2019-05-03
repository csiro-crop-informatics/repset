#!/usr/bin/env bash

biokanga align \
  --sfx ${idxmeta.target}.sfx \
  --format 5 \
  --in ${r1} \
  --pair ${r2}  \
  --out sam \
  --threads ${task.cpus} \
  ${ALIGN_PARAMS}