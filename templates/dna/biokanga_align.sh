#!/usr/bin/env bash

# [species:A_thaliana, version:TAIR10_Pt, simulator:MasonIllumina, nreads:200, mode:PE, length:100, dist:300, distanceDev:50]
# [tool:biokanga, target:A_thaliana_TAIR10_Pt.fasta, seqtype:DNA]


biokanga align \
  --sfx ${idxmeta.target}.sfx \
  --mode 0 \
  --format 5 \
  --pemode 3 \
  --in 1.fq.gz \
  --pair 2.fq.gz  \
  --out out.bam \
  --threads ${task.cpus} \
  --substitutions 5