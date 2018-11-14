# biokanga-manuscript
A repository for a manuscript (application note?)  about [biokanga](https://github.com/csiro-crop-informatics/biokanga).

# Experiments
Adapting and containerising  earlier experiments (ongoing).

On our cluster, running the current version requires close to 60 CPU-days. For a quick test run use `--debug` flag.
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

## Recorded execution

All experiments reported in the manuscript were carried out on a SLURM cluster using:

* Java
```
openjdk version "1.8.0_171"
OpenJDK Runtime Environment (IcedTea 3.8.0) (build 1.8.0_171-b11 suse-27.19.1-x86_64)
OpenJDK 64-Bit Server VM (build 25.171-b11, mixed mode)
```
* Singularity version 2.5.0 -> 2.6.0
* Nextflow version 18.10.1.5003

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

