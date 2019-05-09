#!/usr/bin/env bash

bbmap.sh \
  in=${r1} \
  in2=${r2} \
  threads=${task.cpus} \
  keepnames=t \
  sam=1.3 \
  out=sam \
  ${ALIGN_PARAMS}