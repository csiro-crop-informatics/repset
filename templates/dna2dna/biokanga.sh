#!/usr/bin/env bash

# [species:A_thaliana, version:TAIR10_Pt, simulator:MasonIllumina, nreads:200, mode:PE, length:100, dist:300, distanceDev:50]
# [tool:biokanga, target:A_thaliana_TAIR10_Pt.fasta, seqtype:DNA]


biokanga align \
  --sfx ${ref}.sfx \
  --in ${reads[0]} \
  --pair ${reads[1]}  \
  --out out.bam \
  --threads ${task.cpus} \
  ${ALIGN_PARAMS}