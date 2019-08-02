#!/usr/bin/env bash

GSNAP=\$(awk '{tot+=\$2};END{print tot < 2^32 ? "gsnap" : "gsnapl" }' ${fai});
\${GSNAP} \
  -D genomeDir \
  -d GENOME \
  -A sam \
  --nthreads ${task.cpus} \
  ${ALIGN_PARAMS} \
  <(zcat ${reads[0]}) \
  <(zcat ${reads[1]}) \
  > out.sam