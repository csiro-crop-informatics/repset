  #!/usr/bin/env Rscript

  #args <- commandArgs(TRUE)
  library(reshape2)
  library(ggplot2)

  res<-read.delim('/dev/stdin');
  res2 <- melt(res, id.vars = c("tool", "dist", "distanceDev", "mode", "nreads", "simulator", "species", "version","length"))
  pdf(file="simulatedDNA_summaries.pdf", width=16, height=9);
   ggplot(res2, aes(x=tool, y=value,fill=variable)) +
   geom_bar(stat="identity",position = position_stack(reverse = TRUE)) +
   coord_flip() +
   theme(legend.position = "top") +
   facet_grid(simulator~mode~species);
  dev.off();
