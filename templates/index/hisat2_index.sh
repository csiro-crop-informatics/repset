#!/usr/bin/env bash

hisat2-build --large-index ${ref} ${ref} -p ${task.cpus}