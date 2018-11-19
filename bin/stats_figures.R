#!/usr/bin/env r

library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
#library(hrbrthemes)
library(ggrepel)

stats <- read_csv('/dev/stdin') %>%
 mutate(adapters = case_when(adapters ~ "With adapters",  !adapters  ~ "No adapters")) %>%
 mutate(uniqed = case_when(!uniqed ~ "With secondary/supplementary",  uniqed  ~ "No secondary/supplementary"))

#RUN-TIMES
ggplot(stats)+ aes(tool, aligntime*10^-3/60) +
 geom_point(aes(colour = dataset)) +
 labs(title = "Alignment run times ",
     subtitle = "using 10 logical cores",
     x = "Tool",
     y = "Run time (minutes)") +
 facet_wrap(~adapters)
ggsave(file="runtimes.pdf", width=16, height=9);

ggplot(stats %>% filter(var == "total_read_accuracy", paired == "pairs", uniqed == "With secondary/supplementary"))+
aes(value_dbl, aligntime*10^-3/60,colour=tool) +
geom_point() +
geom_label_repel(aes(label=tool)) +
labs(title = "Accuracy vs alignment run times ",
    subtitle = "using 10 logical cores",
    x = "Accuracy",
    y = "Run time (minutes)") +
guides(label=FALSE, color=FALSE) +
facet_wrap(adapters~dataset)
ggsave(file="accuracy-runtime.pdf", width=16, height=9);

#ALIGNMENT RATES
#Get those that report percentages
stat_perc <- stats %>%
 filter(perc)

cPalette <- c("#000000", "#D55E00", "#999999",  "#009E73")

stats_prop <- stats %>%
 filter(var %in% c("total_read_accuracy",
                   "perc_reads_incorrect",
                   "perc_reads_unaligned",
                   "perc_reads_ambiguous"),
       type == "unique")

ggplot(stats_prop, aes(tool, value_dbl)) +
 geom_bar(aes(fill = var), stat = "identity") +
 facet_wrap(adapters~dataset~uniqed) +
 scale_fill_manual(values=cPalette) +
 labs(title = "Read alignment statistics ",
     subtitle = "uniquely aligned reads",
     x = "Tool",
     y = "Percentage",
     fill = "Alignment classification")
 ggsave(file="align-rates-reads.pdf", width=16, height=9);

stats_prop_b <- stats %>%
 filter(var %in% c("total_bases_accuracy",
                   "perc_bases_incorrect",
                   "perc_bases_unaligned",
                   "perc_bases_ambiguous"),
       type == "unique")


ggplot(stats_prop_b, aes(tool, value_dbl)) +
 geom_bar(aes(fill = var), stat = "identity") +
 facet_wrap(adapters~dataset~uniqed) +
 scale_fill_manual(values=cPalette) +
 labs(title = "Base alignment statistics ",
     subtitle = "uniquely aligned bases",
     x = "Tool",
     y = "Percentage",
     fill = "Alignment classification")
 ggsave(file="align-rates-bases.pdf", width=16, height=9);