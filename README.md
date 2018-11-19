

# Experiments

On our cluster, running pipeline [version 0.5](https://github.com/csiro-crop-informatics/biokanga-manuscript/tree/v0.5)  consumed 56 CPU-days.
See execution [report](https://csiro-crop-informatics.github.io/biokanga-manuscript/report.html)
and [timeline](https://csiro-crop-informatics.github.io/biokanga-manuscript/timeline.html).
This run included each of the input datasets in three replicates. Given the experimental context,
replication does not appear to contribute much, so it may suffice to execute the piepeline with a single replicate using `--replicates 1`.

For a quick test run use `--debug` flag.
In this case only simulated reads from a single dataset and coming from a single human chromosome are aligned to it.
Specific chromosome can be defined using `--debugChromosome ` which defaults to `chr21`.
Additional flag `--adapters` will make a debug run a bit longer but the output results should be slightly more interesting.

## Quick test run

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile docker --debug
```

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile singularity --debug
```

## Full pipeline run

There are a few execution options, all require Nextflow and either Docker or Singularity.
See [nextflow.config](nextflow.config#L22-L47) for available execution profiles, e.g. for local execution this could be


```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile docker
```

or on a SLURM cluster

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile slurm,singularity,singularitymodule
```

 Note that `singularitymodule` profile is used to ensure singularity is available on each execution node by loading an appropriate module. This may not be applicable on your system.

## Experimental pipeline overview

![figures/dag.png](figures/dag.png)

[Same DAG before generalisation of indexing and alignment processes to work with multiple tools](figures/dag-old-colmplex.png)

## Execution environment

All experiments reported in the manuscript were carried out on a SLURM cluster using:

* Java
```
openjdk version "1.8.0_171"
OpenJDK Runtime Environment (IcedTea 3.8.0) (build 1.8.0_171-b11 suse-27.19.1-x86_64)
OpenJDK 64-Bit Server VM (build 25.171-b11, mixed mode)
```
* Singularity version 2.5.0 -> 2.6.0
* Nextflow version 18.10.1.5003

# Adding another aligner


After you have cloned this repository:

1. Add an indexing template to [`templates/`](templates/) subdirectory.
2. Add an alignment template to [`templates/`](templates/) subdirectory.
2. Update [conf/containers.config](conf/containers.config) by specifying a docker hub repository from which an image will be pulled by the pipeline.


## Example

Let's be more specific and follow an example. We will add [bwa](https://github.com/lh3/bwa) - this is just an example

### Add indexing template

```
echo \
'#!/usr/bin/env bash

bwa index -a bwtsw -b 1000000000 ${ref}'  \
> templates/bwa_index.sh
```

Note that `${ref}` is a nextflow variable which will be replaced with the reference FATSA path/filename.


### Add alignment template

```
echo -e \
'#!/usr/bin/env bash

bwa mem -t ${task.cpus} -L 1 ${idxmeta.target} ${r1} ${r2}' \
> templates/bwa_align.sh
```

Applicable nextflow variables resolve as follows :

* `${task.cpus}` - number of logical cpus available to the alignment process
* `${idxmeta.target}` - basename of the index file
* `${r1}` and `${r2}` - path/filenames of paired-end reads

In addition we have lowered the clipping penalty `-L` to increase alignment rates for reads spanning introns.

For many an aligner this would be it, but as it turns out bwa demands matching read names in a given pair, which is incompatible with the benchamrking framework.

We circumvent that by updating the align template [`templates/bwa_align.sh`](templates/bwa_align.sh), such that the name of the second read in pair is adjusted to match the name of the first read, on the fly, just prior to alignment and the changed read names a reverted to their originals immediately after the alignment.

```
bwa mem \
  -t ${task.cpus} \
  -L 1 \
  ${idxmeta.target} \
  ${r1} <(sed 's/b\$/a/' ${r2}) \
  | gawk -vOFS="\\t" '\$1 !~ /^@/ && and(\$2,128) {sub(/a\$/,"b",\$1)};{print}' \
  > sam
```




### Specify container

Insert this code

```
withLabel: bwa {
  container = 'genomicpariscentre/bwa:v0.7.15'
}
```

within the `process { }` block in [conf/containers.config](conf/containers.config).

We opt for docker containers which can also be executed using singularity.
Container images are pulled from docker hub, but nextflow is able to access other registries as well as local images, see relevant nextflow [documentation](https://www.nextflow.io/docs/latest/singularity.html#singularity-docker-hub)





# WRiting:

## Source

Application note is drafted in [RMarkdown](https://rmarkdown.rstudio.com/) in [`biokanga-manuscript.Rmd`](biokanga-manuscript.Rmd) file. RMarkdown is well intergrated in RStudio, but if you'd rather write/edit in a text editor of your choice, here is all that should be required to render the manuscript.

## Bibliography

Among the [alternatives available](https://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html#specifying_a_bibliography) we provisionally opt for BibTeX, see [`references.bib`](references.bib).

## Rendering dependencies

* `R` e.g. on ubuntu `sudo apt apt install r-base-core`
* `pandoc` e.g. on ubuntu `sudo apt install pandoc pandoc-citeproc`
* `LaTeX` e.g. on ubuntu `sudo apt install texlive texlive-latex-extra`
* additional R packages installed and loaded by [`render.R`](render.R)


## Rendering

```
./render.R
```

# Reproductivity of the results

**TODO**

All results presented in the manuscript can be reproduced by executing `nextflow run csiro-crop-informatics/biokanaga-manuscript`.

