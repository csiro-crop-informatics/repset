#!/usr/bin/env r

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

stats <- read_csv('/dev/stdin')
head(stats)

ggplot(stats) +
    aes(aligned, aligntime*10^-3/60,colour=tool) +
    geom_point() +
    geom_label_repel(aes(label=tool)) +
    labs(title = "Accuracy vs alignment run times ",
        subtitle = "using 10 logical cores",
        x = "Aligned",
        y = "Run time (minutes)") +
    guides(label=FALSE, color=FALSE) #+
    #facet_wrap(adapters~dataset)
ggsave(file="realRNA_aligned-runtime.pdf", width=16, height=9);

ggplot(stats) +
    aes(paired, aligntime*10^-3/60,colour=tool) +
    geom_point() +
    geom_label_repel(aes(label=tool)) +
    labs(title = "Accuracy vs alignment run times ",
        subtitle = "using 10 logical cores",
        x = "Aligned as pairs",
        y = "Run time (minutes)") +
    guides(label=FALSE, color=FALSE)
ggsave(file="realRNA_aligned-paired-runtime.pdf", width=16, height=9)