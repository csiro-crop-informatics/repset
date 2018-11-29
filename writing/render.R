#!/usr/bin/env Rscript


location <- "~/local/R_libs/"; dir.create(location, recursive = TRUE  )
if(!require(rmarkdown)){
    install.packages("rmarkdown", lib=location, repos='https://cran.csiro.au')
    library(rmarkdown, lib.loc=location)
}
if(!require(rticles)) {
  install.packages("rticles", type = "source", lib=location, repos='https://cran.csiro.au')
  library(rticles, lib.loc=location)
}
if(!require(bookdown)) {
    install.packages("bookdown", lib=location, repos='https://cran.csiro.au')
    library(bookdown, lib.loc=location)
}
#if(!require(kableExtra)){ #for tables
#  install.packages("kableExtra")
#  library(kableExtra)
#}
#if(!require(revealjs)){ #for slides
#    location <- "~/local/R_libs/"; dir.create(location, recursive = TRUE)
#    install.packages("revealjs", lib=location, repos='https://cran.csiro.au')
#    library(revealjs, lib.loc=location)
#}
rmarkdown::render("biokanga-manuscript.Rmd")

