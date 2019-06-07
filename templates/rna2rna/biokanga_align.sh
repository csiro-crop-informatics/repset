#!/usr/bin/env bash

# [species:A_thaliana, version:TAIR10_Pt, simulator:MasonIllumina, nreads:200, mode:PE, length:100, dist:300, distanceDev:50]
# [tool:biokanga, target:A_thaliana_TAIR10_Pt.fasta, seqtype:DNA]


biokanga align \
  --sfx ${idxmeta.target}.sfx \
  --in 1.fq.gz \
  --pair 2.fq.gz  \
  --out out.sam \
  --threads ${task.cpus} \
  --rptsamseqsthres 100000000 \
  ${ALIGN_PARAMS}


  # - should