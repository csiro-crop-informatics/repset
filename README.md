[![Latest GitHub release](https://img.shields.io/github/release/csiro-crop-informatics/biokanga-manuscript.svg?style=flat-square&logo=github&label=latest%20release)](https://github.com/csiro-crop-informatics/biokanga-manuscript/releases)
[![GitHub commits since latest release](https://img.shields.io/github/commits-since/csiro-crop-informatics/biokanga-manuscript/latest.svg?style=flat-square&logo=github)](https://github.com/csiro-crop-informatics/biokanga-manuscript/releases)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.0-orange.svg)](https://www.nextflow.io/)


# Table of Contents <!-- omit in toc -->
- [Dependencies](#dependencies)
- [Preliminaries](#preliminaries)
  - [Note on terminology (mapping vs alignment)](#note-on-terminology-mapping-vs-alignment)
- [Running the pipeline](#running-the-pipeline)
  - [Execution profiles](#execution-profiles)
    - [Running with docker](#running-with-docker)
    - [Running with singularity](#running-with-singularity)
    - [Running on a SLURM cluster](#running-on-a-slurm-cluster)
    - [Running on AWS batch](#running-on-aws-batch)
  - [Mapping modes](#mapping-modes)
  - [Evaluated mappers](#evaluated-mappers)
- [Capturing results and run metadata](#capturing-results-and-run-metadata)
- [Experimental pipeline overview](#experimental-pipeline-overview)
- [Execution environment](#execution-environment)
- [Adding another mapper](#adding-another-mapper)
    - [Template variables](#template-variables)
      - [Indexing](#indexing)
      - [Mapping](#mapping)
    - [Separate template files (optional)](#separate-template-files-optional)
  - [Non-core mapping parameters (optional)](#non-core-mapping-parameters-optional)
  - [Notes on container specification in `conf/mappers.config`](#notes-on-container-specification-in-confmappersconfig)
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

# Preliminaries

The pipeline consists of several, partly dependent paths
which facilitate the evaluation of mappers using
either DNA-  or RNA-Seq data, either ~~real~~ (temporarily unavailable) or simulated.
The paths can be executed separately or in a single run.
When running separately or re-running the pipeline
the `-resume` flag ensures that previously computed
results (or partial results) are re-used.

Default execution will simulate, align and evaluate reads from a small dataset (a single chromosome from the genome assembly of *A thaliana*).
You can use the `--debug` flag to reduce compute requirements in trial runs on the default data set.

## Note on terminology (mapping vs alignment)

Terms related to read mapping/alignment (and related such as pseudoalignment and quasi mapping) are not used consistently in bioinformatics in a way which would make many a mathematician cringe.
We hereby attempt to strictly follow the common convention by consistently propagating these inconsistencies.
For a much(!) more coherent summary refer to [Lior Pachter's blog post](https://liorpachter.wordpress.com/2015/11/01/what-is-a-read-mapping/).

# Running the pipeline

## Execution profiles

There are several ways to execute the pipeline, each requires Nextflow and either Docker or Singularity.
See [nextflow.config](nextflow.config#L74-L110) for available execution profiles, e.g. for local execution this could be


### Running with docker

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile docker
```

### Running with singularity

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile singularity
```

### Running on a SLURM cluster

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile slurm,singularity,singularitymodule
```

Note:
1. `singularitymodule` profile is used to ensure singularity is available on each execution node by loading an appropriate module.
This may have to be adapted for your system in [nextflow.config](https://github.com/csiro-crop-informatics/biokanga-manuscript/blob/6d0be3f77603f67d13ceeee18de83c517e39db96/nextflow.config#L82).
2. Singularity must also be available on the node where you execute the pipeline, e.g. by running `module load singularity/3.2.1` prior to running the pipeline.


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


## Mapping modes

The workflow incorporates three ways of mapping and evaluating reads, `dna2dna`, `rna2rna`, `rna2dna`
and by default all are executed. To restrict execution to one or two of those,
you can run the workflow with e.g. `--mapmode dna2dna` or `--mapmode rna2rna|rna2dna`.

## Evaluated mappers

An alignment/mapping tool is included in the evaluation if appropriate templates are included
as specified below in [Adding another mapper](#adding-another-mapper).
To execute the workflow for only a subset of the available tools, you can specify e.g.

* `--mappers star` - only evaluate a single tool
* `--mappers 'bwa|bowtie2|biokanga'` - evaluate a subset of tools
* `--mappers '^((?!bwa).)*$'` - evaluate all but this tool

Other regular expressions can be specified to taylor the list of evaluated tools.


# Capturing results and run metadata

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

# Experimental pipeline overview

<!-- ![figures/dag.png](figures/dag.png) -->
<!-- <img src="figures/dag.png" alt="drawing" width="400"/> -->

<!-- For comparison, here is [an earlier version of this graph](figures/dag-old-colmplex.png) -  before indexing and alignment processes were generalised to work with multiple tools. This earlier workflow also excludes evaluation based on real RNA-Seq data. -->

# Execution environment

Execution environment is captured in `runmeta.json`.

# Adding another mapper

A mapper may be included for any or all of the mapping modes (`dna2dna, rna2dna, rna2rna`).
In each case the same indexing template will be used.

After you have cloned this repository **add another entry in [`conf/mappers.config`](conf/mappers.config)**, under

```
params {
  mappersDefinitions = [
    //insert here
  ]
}
```

For example, to add a hypothetical `my_mapper` version `1.0` you might add the following:

```groovy
[
  tool: 'my_mapper',
  version: '1.0',
  container: 'path/to/docker/repository/my_mapper:1.0',
  index: 'my_mapper index --input-fasta ${ref} --output-index ${ref}.idx',
  dna2dna:
  '''
  my_mapper align --index ${ref} \
  -1 ${reads[0]} -2 ${reads[1]} \
  --threads ${task.cpus} \
  ${ALIGN_PARAMS} > out.sam
  '''
],
```

Additional script templates can be added for `rna2rna` and `rna2dna` mapping modes.
Script templates must be wrapped in either single (`'script'`) or triple single (`'''script'''`) quotes.
If you would rather keep the templates in separate files follow [these instructions](#separate-template-files-optional).

### Template variables

Applicable **nextflow** (not bash!) variables resolve as follows:

#### Indexing

* `${task.cpus}` - number of cpu threads available to the indexing process
* `${ref}` - the reference FASTA filename - we use it both to specify the input file and the basename of the generated index

#### Mapping

* `${task.cpus}` - number of logical cpus available to the alignment process
* `${ref}` - basename of the index file (sufficient if aligner uses basename to find multi-file index, otherwise appropriate extension may need to be appended, e.g. `${ref}.idx`).
* `${reads[0]}` and `${reads[1]}` - filenames of paired-end reads
* `${ALIGN_PARAMS}` any additional params passed to the aligner
  * Empty by default but one or more sets of params can be defined in [conf/mapping_params.config](conf/mapping_params.config). When multiple sets of params are specified each set is used in separate execution.


### Separate template files (optional)

If you would rather keep the templates in separate files rather than embedded in [`conf/mappers.config`](conf/mappers.config)
you can place each file under the appropriate template directory:

  * [`templates/index`](templates/index)
  * [`templates/dna2dna`](templates/dna2dna)
  * [`templates/rna2dna`](templates/rna2dna)
  * [`templates/rna2rna`](templates/rna2rna)

and instead of including the script template string directly in [`conf/mappers.config`](conf/mappers.config) as we did above, set

* `rna2dna: true,` which will be resolved to `templates/rna2dna/my_mapper.sh`

or, when using a different file name,

* `rna2dna: 'foo_bar.sh',` which will be resolved to `templates/rna2dna/foo_bar.sh`

See the header of [`conf/mappers.config`](conf/mappers.config) for more details and limitations.



## Non-core mapping parameters (optional)

Add one or more sets of mapping parameters to [conf/mapping_params.config](conf/mapping_params.config),
this is meant for parameter space exploration and should include any fine tuning params while the template
should only include core params essential to mapper execution.



## Notes on container specification in `conf/mappers.config`

You can upload a relevant container image to a docker registry (such as Docker Hub) or locate an existing one [e.g. among  quay biocontainers](https://quay.io/organization/biocontainers). If you opt for an existing one, chose one with a specific version tag and a Dockerfile.
Alternatively, follow our procedure below for [defining per-tool container images and docker automated builds](#per-tool-container-images-and-docker-automated-builds)


We opt for Docker containers which can also be executed using Singularity.
Container images are pulled from Docker Hub, but nextflow is able to access other registries and also local images, see relevant [nextflow documentation](https://www.nextflow.io/docs/latest/singularity.html#singularity-docker-hub)

# Per-tool container images and docker automated builds

Dockerfiles for individual tools used can be found under `dockerfiles/`.
This includes various mappers but also other tools used by the pipeline.
For each tool (or tool-set) we created a docker hub/cloud repository and configured automated builds.

## Setting-up an automated build

Builds can be triggered from branches ~~and tags~~.

The following approach relies on creating a branch for a specific version of a tool.
~~The same can be achieved by simply tagging the relevant commit, but this may result in proliferation of tags while branches can be merged into master and deleted while preserving the history.~~
~~If you'd rather use tags, in (2) change the 'Source type' below to 'Tag' and later tag an appropriate commit using `docker/tool/version` pattern rather than committing to a dedicated branch.~~ (tags can be problematic - if tag is based on version of a tool and container needs to be updated, tags may have to be removed/re-added)

1. Create [Docker Cloud](https://cloud.docker.com/) repo for your tool - *do not link to specific GitHub repo or configure automated build at this stage*, but only *after* it has been created - otherwise the tags for containers built later may be malformed.
2. Link the created a [Docker Cloud](https://cloud.docker.com/) repo with this GitHub repo (go to Builds -> Configure Automated Builds)
3. Add an automated build rule (replace `tool` with the name of the tool).

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

In case the automated build is not triggered for a newly created Docker repo, it may help to delete the Docker repo and repeat steps 1-3 above. Then push some innocuous change to the branch to trigger the build.

If everything works as intended, you may update [conf/containers.config](conf/containers.config) to the new tool version.


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

