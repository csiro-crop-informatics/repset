#!/usr/bin/env Rscript

library(rmarkdown)
library(rticles)
library(bookdown)

rmarkdown::render(Sys.glob("*.Rmd"))