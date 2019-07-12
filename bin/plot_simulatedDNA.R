#!/usr/bin/env r

library(jsonlite)
library(ggplot2)
library(reshape2)

res <- fromJSON("summaries.json", flatten = TRUE)

res2 <- melt(res, na.rm = TRUE, id.vars = c("meta.tool",
                              "meta.dist",
                              "meta.distanceDev",
                              "meta.mode",
                              "meta.coverage",
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
  labs(y = 'number of reads', x = 'tool') +
  geom_bar(stat="identity",position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme(legend.position = "top", legend.title=element_blank()) +
  scale_fill_discrete(breaks=c("results.M_1", "results.M_2", "results.w","results.u"),
                        labels=c("R1 correct", "R2 correct", "Wrong", "Unaligned"))+
  facet_grid(meta.simulator~meta.mode~meta.species)

ggsave(file="simulatedDNA_summaries.pdf", width=16, height=9);