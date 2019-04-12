#!/usr/bin/env r


library(ggplot2)
library(readr)
library(reshape2)

res<-read_csv('/dev/stdin');
res2 <- melt(res, id.vars = c("tool", "dist", "distanceDev", "mode", "nreads", "simulator", "species", "version","length"))
head(res2)
pdf(file="simulatedDNA_summaries.pdf", width=16, height=9);
 ggplot(res2, aes(x=tool, y=value,fill=variable)) +
 geom_bar(stat="identity",position = position_stack(reverse = TRUE)) +
 coord_flip() +
 theme(legend.position = "top") +
 facet_grid(simulator~mode~species);
dev.off();
