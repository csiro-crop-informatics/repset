#!/usr/bin/env bash

bbmap.sh \
  in=${reads[0]} \
  in2=${reads[1]} \
  threads=${task.cpus} \
  keepnames=t \
  sam=1.3 \
  out=out.sam \
  usejni=t \
  Xmx=${task.memory.toMega()}M \
  ${ALIGN_PARAMS}