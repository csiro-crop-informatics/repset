"""
#!/usr/bin/env r

#args <- commandArgs(TRUE)
location <- "~/local/R_libs/"; dir.create(location, recursive = TRUE  )
if(!require(tidyverse)){
  install.packages("tidyverse", lib = location, repos='https://cran.csiro.au')
  library(tidyverse) #, lib.loc = location)
}
#res<-read.delim(gzfile("details.tsv.gz"));
details<-read.delim("details.tsv");

# pdf(file="details.pdf", width=16, height=9);
# binWidth = ${binWidth}
#  details %>%
 #filter(!Chromosome %in% c("Mt","Pt","*","chrUn")) %>%
 #filter(Chromosome %in% c("chr2D")) %>%
 #filter(Class %in% c("w")) %>%
 ggplot(aes(Position, fill=Class)) +
 geom_density(alpha=0.1, bw = ${binWidth}) +
 #geom_vline(xintercept = c(peak), colour="red", linetype="longdash", size=0.5) +
 facet_grid(Species~Aligner~Chromosome~Mode)



  #ggplot(res, aes(x=Position,colour=Class, fill=Class)) +
  #  geom_density(alpha=0.1, adjust=1/10) +
  #  facet_grid(Species~Chromosome~Simulator~Aligner~Mode);
dev.off();

// pdf(file="details7.pdf", width=16, height=9);
details %>%
 # filter(Species %in% c("T_aestivum")) %>% head()
 # filter(!Chromosome %in% c("Mt","Pt","*","chrUn")) %>%
 # filter(str_detect(Chromosome, "^chr1")) %>%
  filter(!Class %in% c("u")) %>%
  filter(Simulator %in% c("MasonIllumina")) %>%
ggplot(aes(x=Position, fill = Class, colour=Class)) +
  geom_density(aes(x=Position, y=..count..*${binWidth}), alpha=0.1, bw = ${binWidth}) +
  facet_grid(Species ~ Chromosome ~ Aligner  ~ Mode)
// dev.off();
"""