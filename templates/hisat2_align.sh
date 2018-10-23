#!/usr/bin/env bash

#TUNING_PARAMS="-N 1 -L 20 -i S,1,0.5 -D 25 -R 5 --pen-noncansplice 12  --mp 1,0  --sp 3,0"

hisat2 -x ${idxmeta.target} \
  -1 ${r1} \
  -2 ${r2} \
  --threads ${task.cpus} \
  --reorder \
  -f \
  > sam