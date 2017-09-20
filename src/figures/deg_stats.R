
library("readr")
library("ggplot2")
library("dplyr")
library("tidyr")
library("ggthemes")

dir.create("figures/deg_stats")

# setwd("../..")
deg_stats <- read_tsv("results/DEG_stats_by_group.tsv")

tit <- c("EarlyDengue vs Control", "LateDengue vs Control")
s <- deg_stats %>% 
    filter(title %in% tit, Direction %in% c("Down-regulated", "Up-regulated")) %>% 
    mutate(n_sig = n*ifelse(Direction=="Down-regulated", -1, 1))
ggplot(s, aes(x=Group, fill=Direction, y=n_sig)) +
       geom_bar(stat="identity") +
       geom_text(aes(label=n, vjust=-1*sign(n_sig)), size=7) +
       scale_fill_manual(values=c("Down-regulated"="darkgreen", "Up-regulated"="darkred")) +
       theme_gdocs() +
       facet_grid(. ~ title)
ggsave("figures/deg_stats/deg_stats_groups.pdf", width=12)
