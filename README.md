[![Latest GitHub tag](https://img.shields.io/github/tag/csiro-crop-informatics/biokanga-manuscript.svg?label=latest%20release&logo=github&style=flat-square)](https://github.com/csiro-crop-informatics/biokanga-manuscript/releases)

# Table of Contents <!-- omit in toc -->
- [Experiments](#experiments)
  - [Quick test run](#quick-test-run)
    - [Running nextflow with singularity](#running-nextflow-with-singularity)
    - [Running nextflow with docker](#running-nextflow-with-docker)
    - [Running on AWS batch](#running-on-aws-batch)
  - [Full pipeline run](#full-pipeline-run)
  - [Experimental pipeline overview](#experimental-pipeline-overview)
  - [Execution environment](#execution-environment)
- [Adding another aligner](#adding-another-aligner)
  - [Example](#example)
    - [Add indexing template](#add-indexing-template)
    - [Add alignment template](#add-alignment-template)
    - [Specify container](#specify-container)
- [WRiting](#writing)
  - [Rendering outside the pipeline](#rendering-outside-the-pipeline)
    - [Using docker](#using-docker)
    - [Using singularity](#using-singularity)
    - [Natively](#natively)
  - [Bibliography](#bibliography)
# Experiments

On our cluster, running pipeline [version 0.5](https://github.com/csiro-crop-informatics/biokanga-manuscript/releases/tag/v0.5) consumed 56 CPU-days.
See execution [report](https://csiro-crop-informatics.github.io/biokanga-manuscript/report.html)
and [timeline](https://csiro-crop-informatics.github.io/biokanga-manuscript/timeline.html).
This run included each of the input datasets in three replicates. Given the experimental context,
replication does not appear to contribute much, so it may suffice to execute the pipeline with a single replicate using `--replicates 1`,
thus reducing the CPU-time to under 8 days (based on a run of [version 0.6](https://github.com/csiro-crop-informatics/biokanga-manuscript/releases/tag/v0.6)).

For a quick test run use `--debug` flag.
In this case only simulated reads from a single dataset and coming from a single human chromosome are aligned to it.
Specific chromosome can be defined using `--debugChromosome ` which defaults to `chr21`. By default, all pre-defined aligners are executed.
To only specify a single aligner you can e.g. use `--aligners biokanga` or for several aligners e.g. `--aligners 'biokanga|dart|hisat2'`.

Additional flag `--adapters` will make a debug run a bit longer but the output results should be slightly more interesting by including datsets with retained adapters.


## Quick test run

### Running nextflow with singularity

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile docker --debug --aligners 'biokanga|hisat2|star'
```

### Running nextflow with docker

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile singularity --debug --aligners 'biokanga|hisat2|star'
```

### Running on AWS batch

If you are new to AWS batch and/or nextflow, follow [this blog post](https://antunderwood.gitlab.io/bioinformant-blog/posts/running_nextflow_on_aws_batch/), once you are done, or you already use AWS batch, simply run

```
nextflow run csiro-crop-informatics/biokanga-manuscript \
  -profile awsbatch --debug --aligners 'biokanga|hisat2|star' \
  -work-dir s3://your_s3_bucket/work --outdir s3://your_s3_bucket/results
```

after replacing `your_s3_bucket` with a bucket you have created on S3.

[**Warning! You will be charged by AWS according to your resource use.**](https://docs.aws.amazon.com/awsaccountbilling/latest/aboutv2/monitoring-costs.html)


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

To run the pipeline on [AWS batch](https://aws.amazon.com/batch/), follow the [instructions above](#running-on-aws-batch) but drop the `--debug` flag.

## Experimental pipeline overview

![figures/dag.png](figures/dag.png)

For comparison, here is [an earlier version of this graph](figures/dag-old-colmplex.png) -  before indexing and alignment processes were generalised to work with multiple tools.

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
3. Update [conf/containers.config](conf/containers.config) by specifying a docker hub repository from which an image will be pulled by the pipeline.


## Example

Let's be more specific and follow an example. We will add [bowtie2](https://github.com/BenLangmead/bowtie2/releases/tag/v2.3.4.1).

### Add indexing template

```
echo \
'#!/usr/bin/env bash

bowtie2-build --threads ${task.cpus} ${ref} ${ref}
> templates/bowtie2_index.sh
```

Applicable nextflow variables resolve as follows:

* `${task.cpus}` - number of cpu threads available to the alignment process
* `${ref}` - the reference FASTA path/filename - in this case we use it both to specify the input file and the basename of the generated index


### Add alignment template

```
echo -e \
'#!/usr/bin/env bash

bowtie2 \
  -p ${task.cpus} \
  -x ${idxmeta.target} \
  -1 ${r1} \
  -2 ${r2} \
  -f \
  --threads  ${task.cpus} \
  --local \
  > sam' \
> templates/bowtie2_align.sh
```

Applicable nextflow variables resolve as follows :

* `${task.cpus}` - number of logical cpus available to the alignment process
* `${idxmeta.target}` - basename of the index file
* `${r1}` and `${r2}` - path/filenames of paired-end reads

In addition we have used bowtie's `--local` flag to increase alignment rates for reads spanning introns.


### Specify container


1. Upload a relevant container image to docker hub or [locate an existing one](https://hub.docker.com/search/?isAutomated=0&isOfficial=0&page=1&pullCount=0&q=bowtie2&starCount=0). If you opt for an existing one, chose one with a specific version tag and a Dockerfile.

2. Insert container specification

```
withLabel: bowtie2 {
  container = 'comics/bowtie2:2.3.4.1'
}
```
within the `process {   }` block in [conf/containers.config](conf/containers.config).

We opt for docker containers which can also be executed using singularity.
Container images are pulled from docker hub, but nextflow is able to access other registries and also local images, see relevant [nextflow documentation](https://www.nextflow.io/docs/latest/singularity.html#singularity-docker-hub)


# WRiting

Application note is drafted in [RMarkdown](https://rmarkdown.rstudio.com/) in [`writing/biokanga-manuscript.Rmd`](writing/biokanga-manuscript.Rmd) file.
RMarkdown is well integrated in RStudio, but can be written/edited in a text editor of your choice.
Rendering of the manuscript constitutes the final step of our nextflow pipeline which relies on a container defined in [`dockerfiles/renderer.Dockerfile`](dockerfiles/renderer.Dockerfile) for rendering environment.

## Rendering outside the pipeline

### Using docker

```sh
docker run --rm --user $(id -u):$(id -g) \
  --volume $(pwd)/writing:/writing \
  --workdir /writing rsuchecki/renderer:0.1 ./render.R
```

### Using singularity

```sh
singularity exec --pwd $(pwd)/writing docker://rsuchecki/renderer:0.1 ./render.R
```

### Natively

If you'd like to render the manuscript without docker/singularity, you will need the following:

* `R` e.g. on ubuntu `sudo apt apt install r-base-core`
* `pandoc` e.g. on ubuntu `sudo apt install pandoc pandoc-citeproc`
* `LaTeX` e.g. on ubuntu `sudo apt install texlive texlive-latex-extra`
* `R` packages:
  * `rmarkdown`
  * `rticles`
  * `bookdown`

Then:

```
cd writing && ./render.R
```

## Bibliography

Among the [alternatives available](https://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html#specifying_a_bibliography) we opted for BibTeX, see [`writing/references.bib`](writing/references.bib).


