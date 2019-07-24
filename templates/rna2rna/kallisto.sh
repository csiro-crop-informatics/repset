#!/usr/bin/env bash

kallisto quant \
  --threads ${task.cpus} \
  ${ALIGN_PARAMS} \
  --index ${ref}.idx \
  --pseudobam \
  --output-dir . \
  ${reads}
#&& mv output/pseudoalignments.bam out.bam

