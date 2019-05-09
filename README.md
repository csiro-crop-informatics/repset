[![Latest GitHub release](https://img.shields.io/github/release/csiro-crop-informatics/biokanga-manuscript.svg?style=flat-square&logo=github&label=latest%20release)](https://github.com/csiro-crop-informatics/biokanga-manuscript/releases)
[![GitHub commits since latest release](https://img.shields.io/github/commits-since/csiro-crop-informatics/biokanga-manuscript/latest.svg?style=flat-square&logo=github)](https://github.com/csiro-crop-informatics/biokanga-manuscript/releases)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.0-orange.svg)](https://www.nextflow.io/)


# Table of Contents <!-- omit in toc -->
- [Dependencies](#dependencies)
- [Experiments](#experiments)
  - [Quick start](#quick-start)
  - [Simulated RNA-Seq](#simulated-rna-seq)
    - [Quick test run](#quick-test-run)
  - [Real RNA-Seq](#real-rna-seq)
  - [Simulated DNA-Seq](#simulated-dna-seq)
  - [Real DNA-Seq](#real-dna-seq)
  - [Execution environments](#execution-environments)
    - [Running nextflow with docker](#running-nextflow-with-docker)
    - [Running nextflow with singularity](#running-nextflow-with-singularity)
    - [Running on AWS batch](#running-on-aws-batch)
  - [Full pipeline run](#full-pipeline-run)
    - [Capturing results and run metadata](#capturing-results-and-run-metadata)
  - [Experimental pipeline overview](#experimental-pipeline-overview)
  - [Execution environment](#execution-environment)
- [Adding another aligner](#adding-another-aligner)
  - [Example](#example)
    - [Add indexing template](#add-indexing-template)
    - [Add RNA alignment template](#add-rna-alignment-template)
    - [Add DNA alignment template](#add-dna-alignment-template)
    - [Add additional alignment parameters](#add-additional-alignment-parameters)
    - [Specify container](#specify-container)
- [Per-tool container images and docker automated builds](#per-tool-container-images-and-docker-automated-builds)
  - [Setting-up an automated build](#setting-up-an-automated-build)
  - [Adding or updating a Dockerfile](#adding-or-updating-a-dockerfile)
- [Report](#report)
  - [Rendering outside the pipeline](#rendering-outside-the-pipeline)
    - [Using docker](#using-docker)
    - [Using singularity](#using-singularity)
    - [Natively](#natively)
- [Manuscript](#manuscript)
  - [Rendering](#rendering)
  - [Bibliography](#bibliography)

# Dependencies

* Nextflow [![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.0-orange.svg)](https://www.nextflow.io/)
* and
  * either Singularity [![Singularity](https://img.shields.io/badge/Singularity-%E2%89%A53.1.1-orange.svg)](https://www.sylabs.io/singularity/)
   * or Docker

# Experiments

The pipeline consists of several, partly dependent paths
which facilitate the evaluation of aligners using
either DNA-  or RNA-Seq data, either real or simulated.
The paths can be executed separately or in a single run.
When running separately or re-running the pipeline
the `-resume` flag ensures that previously computed
results (or partial results) are reused.

## Quick start

Use the `--debug` flag to run the whole pipeline or part of it with reduced input data.

There are several ways to execute the pipeline, each requires Nextflow and either Docker or Singularity.

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile docker --debug
```

or on a SLURM cluster

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile slurm,singularity,singularitymodule --debug
```

Note:
1. `singularitymodule` profile is used to ensure singularity is available on each execution node by loading an appropriate module.
This may need to be adapted for your system in [nextflow.config](https://github.com/csiro-crop-informatics/biokanga-manuscript/blob/6d0be3f77603f67d13ceeee18de83c517e39db96/nextflow.config#L82).
2. Singularity must also be available on the node where you execute the pipeline, e.g. by running `module load singularity/3.1.1` prior to running the pipeline.
3.



## Simulated RNA-Seq
On our cluster, running pipeline [version 0.5](https://github.com/csiro-crop-informatics/biokanga-manuscript/releases/tag/v0.5) consumed 56 CPU-days.
See execution [report](https://csiro-crop-informatics.github.io/biokanga-manuscript/report.html)
and [timeline](https://csiro-crop-informatics.github.io/biokanga-manuscript/timeline.html).
This run included each of the input datasets in three replicates. Given the experimental context,
replication does not appear to contribute much, so it may suffice to execute the pipeline with a single replicate using `--replicates 1`,
thus reducing the CPU-time to under 8 days (based on a run of [version 0.6](https://github.com/csiro-crop-informatics/biokanga-manuscript/releases/tag/v0.6)).

### Quick test run

For a quick test run use the `--debug` flag.
In this case only simulated reads from a single dataset and coming from a single human chromosome are aligned to it.
Specific chromosome can be defined using `--debugChromosome ` which defaults to `chr21`. By default, all pre-defined aligners are executed.
To only specify a single aligner you can e.g. use `--aligners biokanga` or for several aligners e.g. `--alignersRNA 'biokanga|dart|hisat2'`.

Additional flag `--adapters` will make a debug run a bit longer but the output results should be slightly more interesting by including datasets with retained adapters.



## Real RNA-Seq

TODO

## Simulated DNA-Seq

TODO


## Real DNA-Seq

TODO


## Execution environments

We provide several execution profiles, as before, you may use the `--debug` to test the pipeline and the execution environment before a full run.

### Running nextflow with docker

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile docker
```

### Running nextflow with singularity

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile singularity
```

### Running on AWS batch

If you are new to AWS batch and/or nextflow, follow [this blog post](https://antunderwood.gitlab.io/bioinformant-blog/posts/running_nextflow_on_aws_batch/), once you are done, or you already use AWS batch, simply run

```
nextflow run csiro-crop-informatics/biokanga-manuscript \
  -profile awsbatch \
  -work-dir s3://your_s3_bucket/work \
  --outdir s3://your_s3_bucket/results
```

after replacing `your_s3_bucket` with a bucket you have created on S3.

[**Warning! You will be charged by AWS according to your resource use.**](https://docs.aws.amazon.com/awsaccountbilling/latest/aboutv2/monitoring-costs.html)


## Full pipeline run

There are a few ways to execute the pipeline, all require Nextflow and either Docker or Singularity.
See [nextflow.config](nextflow.config#L46-L84) for available execution profiles, e.g. for local execution this could be

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile docker
```

or on a SLURM cluster

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile slurm,singularity,singularitymodule
```

Note that `singularitymodule` profile is used to ensure singularity is available on each execution node by loading an appropriate module.
This may need to be adapted for your system.
In addition Singularity must also be available on the node where you execute the pipeline.

To run the pipeline on [AWS batch](https://aws.amazon.com/batch/), follow the [instructions above](#running-on-aws-batch).

### Capturing results and run metadata

Each pipeline run generates a number of files including
* results in the form of report, figures, tables etc.
* run metadata reflecting information about the pipeline version, software and compute environment etc.

These can be simply collected from the output directories but for full traceability of the results, the following procedure is preferable:
1. Select a tagged revision or add a tag (adhering to the usual semantic versioning approach)
2. Generate a [Git Hub access token](https://help.github.com/en/articles/creating-a-personal-access-token-for-the-command-line)
   which will allow the pipeline to create releases in this or a forked repository,
   when creating the token it suffices to select only the following scope:
   > `public_repo`   Access public repositories
3. Make the access token accessible as an environmental variable e.g. `GH_TOKEN='your-token-goes-here'`
4. Run the pipeline from the remote repository, specifying
    - the required revision  e.g. `-revision v0.8.3`
    - the `--release` flag

On successful completion of the pipeline a series of API calls will be made to

1. create a new release
2. upload results and metadata files as artefacts for that release
3. finalize the release

The last of this calls will trigger minting of a DOI for that release if Zenodo integration is configured and enabled for the repository.

## Experimental pipeline overview

![figures/dag.png](figures/dag.png)

<!-- For comparison, here is [an earlier version of this graph](figures/dag-old-colmplex.png) -  before indexing and alignment processes were generalised to work with multiple tools. This earlier workflow also excludes evaluation based on real RNA-Seq data. -->

## Execution environment

Execution environment is captured in `runmeta.json`.

# Adding another aligner

An aligner may be included for DNA alignment, RNA alignment or both. In each case the same indexing template will be used.

After you have cloned this repository:

1. Add an indexing template to [`templates/index`](templates/index) subdirectory.
2. Add an alignment template(s) to [`templates/rna`](templates/rna) and/or [`templates/dna`](templates/dna) subdirectories.
3. Add one or more sets of aligner parameters to [conf/aligners.config](conf/aligners.config).
4. Update [conf/containers.config](conf/containers.config) by specifying a docker hub repository from which an image will be pulled by the pipeline.


## Example

Let's be more specific and follow an example. We will add [bowtie2](https://github.com/BenLangmead/bowtie2/releases/tag/v2.3.5).

### Add indexing template

```
echo \
'#!/usr/bin/env bash

bowtie2-build --threads ${task.cpus} ${ref} ${ref}
> templates/index/bowtie2_index.sh
```

Applicable nextflow variables resolve as follows:

* `${task.cpus}` - number of cpu threads available to the alignment process
* `${ref}` - the reference FASTA path/filename - in this case we use it both to specify the input file and the basename of the generated index


### Add RNA alignment template

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
> templates/rna/bowtie2_align.sh
```

Applicable nextflow variables resolve as follows :

* `${task.cpus}` - number of logical cpus available to the alignment process
* `${idxmeta.target}` - basename of the index file
* `${r1}` and `${r2}` - path/filenames of paired-end reads

In addition we have used bowtie's `--local` flag to increase alignment rates for reads spanning introns.

### Add DNA alignment template

TODO

### Add additional alignment parameters



### Specify container

1. Upload a relevant container image to docker hub or [locate an existing one](https://hub.docker.com/search/?isAutomated=0&isOfficial=0&page=1&pullCount=0&q=bowtie2&starCount=0). If you opt for an existing one, chose one with a specific version tag and a Dockerfile.
Alternatively, follow our procedure below for [defining per-tool container images and docker automated builds](#per-tool-container-images-and-docker-automated-builds)

2. Insert container specification

```
withLabel: bowtie2 {
  quay.io/biocontainers/bowtie2:2.3.5--py27he860b03_0
}
```
within the `process {   }` block in [conf/containers.config](conf/containers.config).

We opt for docker containers which can also be executed using singularity.
Container images are pulled from docker hub, but nextflow is able to access other registries and also local images, see relevant [nextflow documentation](https://www.nextflow.io/docs/latest/singularity.html#singularity-docker-hub)

# Per-tool container images and docker automated builds

Dockerfiles for individual tools used can be found under `dockerfiles/`.
This includes various aligners but also other tools used by the pipeline.
For each tool we created a docker hub/cloud repository and configured automated builds.

## Setting-up an automated build

Builds can be triggered from branches ~~and tags~~.

The following approach relies on creating a branch for a specific version of a tool.
~~The same can be achieved by simply tagging the relevant commit, but this may result in proliferation of tags while branches can be merged into master and deleted while preserving the history.~~
~~If you'd rather use tags, in (2) change the 'Source type' below to 'Tag' and later tag an appropriate commit using `docker/tool/version` pattern rather than committing to a dedicated branch.~~ (tags can be problematic - if tag is based on version of a tool and container needs to be updated, tags may have to be removed/re-added)

1. Link a [Docker Cloud](https://cloud.docker.com/) repo with this GitHub repo (go to Builds -> Configure Automated Builds)
2. Add an automated build rule (replace `tool` with the name of the tool).

| Source type   | Source                   | Docker Tag  | Dockerfile location | Build Context  |
| ------------- | ------------------------ | ----------- | ------------------- | -------------- |
| Branch        | `/^docker\/tool\/(.*)$/` | `{\1}`      | `tool.Dockerfile`   | `/dockerfiles` |


## Adding or updating a Dockerfile

Checkout a new branch replacing `tool` and `version` with the intended tool name and version, respectively.
For example,

```
tool='bwa'
version='0.7.17'
```

```
git checkout -b docker/${tool}/${version}
```

Add or modify `dockerfiles/${tool}.Dockerfile` as required.

Commit and push to trigger an automated build

```
git add dockerfiles/${tool}.Dockerfile
git commit
git push --set-upstream origin docker/${tool}/${version}
```

This should trigger an automated build in the linked Docker Hub/cloud repository.
If everything works as intended, you may update [conf/containers.config](conf/containers.config) to the new tool version

Then either create a PR to merge the new branch into master or,
if you have write permissions for this repository or working on your fork of it, checkout master and merge.

```
git checkout master
git merge docker/${tool}/${version}
```

# Report

TODO: add information on

* how to edit the report template
* how the final report gets generated

If report template is sufficiently generic we will be able to easily render to html and PDF, otherwise we should settle for HTML(?).

Rendering of the report constitutes the final step of the pipeline and relies on a container defined in [`dockerfiles/renderer.Dockerfile`](dockerfiles/renderer.Dockerfile) for rendering environment.

## Rendering outside the pipeline

There are several ways for rendering of the report outside the pipeline, with docker being the preferred option.

### Using docker

```sh
docker run --rm --user $(id -u):$(id -g) \
  --volume $(pwd)/report:/report \
  --workdir /report rsuchecki/renderer:0.2 ./render.R
```

### Using singularity

```sh
singularity exec --pwd $(pwd)/report docker://rsuchecki/renderer:0.1 ./render.R
```

### Natively

If you'd like to render the report without docker/singularity, you will need the following:

* `R` e.g. on ubuntu `sudo apt apt install r-base-core`
* `pandoc` e.g. on ubuntu `sudo apt install pandoc pandoc-citeproc`
* `LaTeX` e.g. on ubuntu `sudo apt install texlive texlive-latex-extra`
* `R` packages:
  * `rmarkdown`
  * `rticles`
  * `bookdown`

Then:

```
cd report && ./render.R
```



# Manuscript

Manuscript source is under `manuscript/` sub directory on `manuscript` branch which should not be merged into master.
Application note is drafted in [RMarkdown](https://rmarkdown.rstudio.com/) in [`manuscript/biokanga-manuscript.Rmd`](blob/manuscript/manuscript/biokanga-manuscript.Rmd) file.
RMarkdown is well integrated in RStudio, but can be written/edited in a text editor of your choice.

## Rendering
The manuscript will be rendered the pipeline is executed while `manuscript` branch is checked out, either
  * locally
  or
  * by specifying `-revision manuscript` at run-time

Appropriate revision of the master branch should first be mnerged into the `manuscript` branch.

The manuscript can be rendered outside the pipeline in a fashion analogous to how this can be done for the [report](#rendering-outside-the-pipeline),
just replace any use of `report` by `manuscript`.





## Bibliography

Among the [alternatives available](https://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html#specifying_a_bibliography) we opted for BibTeX, see [`writing/references.bib`](writing/references.bib).

