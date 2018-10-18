#!/usr/bin/env Rscript

location <- "~/local/R_libs/"; dir.create(location, recursive = TRUE  )
if(!require(tidyverse)){
  install.packages("tidyverse", lib = location, repos='https://cran.csiro.au')
  library(tidyverse, lib.loc = location)
}

new_vars_single <- c("total_reads",
              "total_read_accuracy",
              "unique_read_accuracy",
              "perc_reads_incorrect",
              "perc_reads_ambiguous",
              "perc_reads_unaligned",
              "perc_reads_aligned",
              "perc_reads_true_introns",
              "total_bases",
              "total_bases_accuracy",
              "unique_bases_accuracy",
              "perc_bases_incorrect",
              "perc_bases_ambiguous",
              "perc_bases_unaligned",
              "perc_bases_aligned",
              "perc_bases_insertions",
              "perc_bases_deletions",
              "insertion_fd",
              "insertion_fn",
              "deletion_fd",
              "deletion_fn",
              "skipping_fd",
              "skipping_fn",
              "junction_fd",
              "junction_fn")

new_vars_multi <- c("total_reads",
                    "total_read_accuracy",
                    "unique_read_accuracy",
                    "perc_reads_incorrect",
                    "perc_reads_ambiguous",
                    "perc_reads_unaligned",
                    "perc_reads_aligned",
                    "total_bases",
                    "total_bases_accuracy",
                    "unique_bases_accuracy",
                    "perc_bases_incorrect",
                    "perc_bases_ambiguous",
                    "perc_bases_unaligned",
                    "perc_bases_aligned")

#Handle multimappers col type: single/multi

parse_file <- function(filename, multi = FALSE, debug = FALSE) {
  if(debug) browser()
  if(grepl("multi_mappers", filename)) multi <- TRUE
  raw <- read_tsv(filename, comment = "-", col_names = c("var", "value"))

  #deal with | sep rows (26:27)
  junctions <- raw %>%
    filter(grepl("Junctions", var)) %>%
    separate(value, paste0("junction_", c("none", "left", "right", "both")), "\\|") %>%
    mutate(var = c("n", "perc")) %>%
    gather(var_new, value, -var) %>%
    mutate(var = paste(var_new, var, sep = "_"),
           paired = "single") %>%
    select(-var_new)

  if(multi) new_vars <- c(new_vars_single, new_vars_multi) else new_vars <- new_vars_single
  paired <- rep("single", 25)
  if(multi) paired <- c(paired, rep("pairs", 14))

  other <- raw %>%
    filter(!grepl("Junctions", var)) %>%
    select(-var) %>%
    mutate(var = new_vars,
           paired = paired)

  out <- bind_rows(other, junctions)

  out <- out %>%
    mutate(perc = ifelse(grepl("%", value), TRUE, FALSE),
           value_dbl = ifelse(grepl("no.*called", value), NA,
                                 as.numeric(sub("([^%\\(]*).*", "\\1", value))))
}

write_csv(parse_file('/dev/stdin'),'/dev/stdout')