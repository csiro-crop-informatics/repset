#!/usr/bin/env bash

kallisto quant \
  --threads ${task.cpus} \
  ${ALIGN_PARAMS} \
  --index ${idxmeta.target}.idx \
  --pseudobam \
  --output-dir . \
  ${reads}
#&& mv output/pseudoalignments.bam out.bam

