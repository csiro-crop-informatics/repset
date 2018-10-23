#!/usr/bin/env bash

biokanga align \
  --sfx ${idxmeta.target}.sfx \
  --mode 0 \
  --format 5 \
  --maxns 2 \
  --pemode 3 \
  --pairmaxlen 50000 \
  --in ${r1} \
  --pair ${r2}  \
  --out sam \
  --threads ${task.cpus} \
  --substitutions 5 \
  --minchimeric 50