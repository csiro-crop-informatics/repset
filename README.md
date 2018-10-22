# biokanga-manuscript
A repository for a manuscript (application note?)  about [biokanga](https://github.com/csiro-crop-informatics/biokanga).

# Experiments
Adapting and containerising  earlier experiments (ongoing).

Running the current version requires approximately 220 CPU-hours. For a quick-ish test run use `--debug` flag.
In this case only the simulated reads coming from a single human chromosome are aligned to it.
 Specific chromosome can be defined using `--debugChromosome ` which defaults to `chr21`.

There are a few execution options, all require Nextflow and Singularity.

```
nextflow run csiro-crop-informatics/biokanga-manuscript
```

See [nextflow.config](nextflow.config#L22-L40) for available execution profiles, e.g. for local execution this could be

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile singularity
```

or on a SLURM cluster

```
nextflow run csiro-crop-informatics/biokanga-manuscript -profile slurm,singularity,singularitymodule
```

## Experimental pipeline overview


![figures/dag.png](figures/dag.png)
[Same DAG before generalisation of indexing and alignment processes to work with multiple tools](figures/dag-old-colmplex.png)


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

**TODO: containerize the rendering environment**

## Rendering

```
./render.R
```

# Reproductivity of the results

**TODO**

All results presented in the manuscript can be reproduced by executing `nextflow run csiro-crop-informatics/biokanaga-manuscript`.

