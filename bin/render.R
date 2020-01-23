#!/usr/bin/env Rscript

library(rmarkdown)
library(rticles)
library(bookdown)
library(tidyverse)
library(jsonlite)
library(kableExtra)

rmarkdown::render(Sys.glob("*.Rmd"))