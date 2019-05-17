#!/usr/bin/env r

library(jsonlite)
library(ggplot2)
library(reshape2)

res <- fromJSON("summaries.json", flatten = TRUE)

res2 <- melt(res, na.rm = TRUE, id.vars = c("meta.tool",
                              "meta.dist",
                              "meta.distanceDev",
                              "meta.mode",
                              "meta.nreads",
                              "meta.simulator",
                              "meta.species",
                              "meta.version",
                              "meta.length",
                              "meta.paramslabel",
                              "meta.ALIGN_PARAMS",
                              "meta.seqtype",
                              "meta.target"))

ggplot(na.omit(res2), aes(x=meta.tool, y=as.numeric(value), fill=variable)) +
  scale_y_continuous()+
  geom_bar(stat="identity",position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme(legend.position = "top") +
  facet_grid(meta.simulator~meta.mode~meta.species)

ggsave(file="simulatedDNA_summaries.pdf", width=16, height=9);