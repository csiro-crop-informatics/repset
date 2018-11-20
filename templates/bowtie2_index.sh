#!/usr/bin/env bash

bowtie2-build --threads ${task.cpus} ${ref} ${ref}