#!/usr/bin/env bash

bbmap.sh \
  in=${reads[0]} \
  in2=${reads[1]} \
  threads=${task.cpus} \
  keepnames=t \
  sam=1.3 \
  out=sam \
  usejni=t \
  ${ALIGN_PARAMS}