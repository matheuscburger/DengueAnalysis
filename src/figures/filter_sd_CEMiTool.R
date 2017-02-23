#!/usr/bin/env Rscript

library("tidyverse")
library("stringr")
library("pracma")

dir.create(file.path('figures', 'filter_sd_cemitool'))

expression_files_directory <- "data/processed/collapsed"
reannotation_filename <- "config/reannotation/annotation_long.tsv"

expr_filenames <- list.files(expression_files_directory, full.names=T)
names(expr_filenames) <- expr_filenames %>% str_match(".+/(.+)\\.tsv") %>% .[, 2]

expr <- tibble(Symbol=character(0),
               Sample=character(0),
               Expression=numeric(0),
               Study=character(0))
for(study_name in names(expr_filenames)) {
    filename <- expr_filenames[study_name]
    expr <- read_tsv(filename) %>% 
        gather("Sample", "Expression", -Symbol) %>% 
        mutate(Study=study_name) %>%
        union(expr)
}

sd_genes <- expr %>%
    group_by(Symbol, Study) %>%
    summarise(SD=sd(Expression)) %>%
    group_by(Study) %>%
    arrange(SD) %>%
    mutate(RankSD=length(SD)-rank(SD, ties.method="first")) %>%
    ungroup()

pl <- ggplot(sd_genes, aes(x=SD, y=RankSD, group=Study, color=Study)) + 
    geom_line(size=1) +
    labs(y="#Genes", x="Standard Deviation", title="Gene Filtering") +
    xlim(0, 1.5) +
    theme_minimal()
ggsave(pl, filename=file.path("figures", "filter_sd_cemitool", "ngenes_sd.pdf"))

## by class

reannot <- read_tsv(reannotation_filename) %>%
    select(Gene, Type, Group) %>%
    unique

sd_extended <- sd_genes %>% 
    left_join(reannot, by=c("Symbol"="Gene")) %>%
    group_by(Study, Group) %>%
    arrange(SD) %>%
    mutate(RankSD=length(SD)-rank(SD, ties.method="first")) %>%
    ungroup()

pl <- ggplot(sd_extended, aes(x=SD, y=RankSD, group=Study, color=Study)) + 
    geom_line(size=1) +
    facet_wrap(~Group) +
    labs(y="#Genes", x="Standard Deviation", title="Gene Filtering") +
    xlim(0, 1.5) +
    theme_minimal()
ggsave(pl, filename=file.path("figures", "filter_sd_cemitool", "allclasses_sd.pdf"), width=14, height=14)

pl <- sd_extended %>%
    dplyr::filter(Group=="lnoncoding") %>% 
    ggplot(., aes(x=SD, y=RankSD, group=Study, color=Study)) +
    geom_line(size=1) +
    labs(y="#Genes", x="Standard Deviation", title="Gene Filtering\nLong Noncoding RNA") +
    xlim(0, 1.5) +
    theme_minimal()
ggsave(pl, filename=file.path("figures", "filter_sd_cemitool", "lnoncoding_sd.pdf"))



# number of lncRNAs in two platforms

study2plat <- list("GSE13052" = "GPL2700",
                   "GSE28405" = "GPL2700",
                   "GSE51808"="GPL570",
                   "GSE43777"="GPL570")
lnc_in_plat <- sd_extended %>% mutate(Platform=unlist(study2plat[Study])) %>% 
    dplyr::filter(Group=="lnoncoding") %>%
    dplyr::select(Symbol, Platform, Group, Type) %>% 
    mutate(Present=TRUE) %>% 
    unique %>%
    spread(key=Platform, value=Present, convert=T, fill=FALSE) %>%
    mutate(Both=GPL2700 & GPL570)

sum(lnc_in_plat[["Both"]])

in_both <- lnc_in_plat %>% filter(Both) %>% select(Symbol) %>% .[[1]]

## tentando fazer um gráfico juntando os 4 estudos
sd_min <- sd_extended %>%
    filter(Group=="lnoncoding", Symbol %in% in_both) %>%
    ungroup() %>%
    group_by(Symbol) %>%
    mutate(MinSD=min(SD)) %>%
    ungroup() %>% 
    select(-Study, -SD, -RankSD) %>%
    unique() %>%
    group_by(Group) %>%
    arrange(MinSD) %>%
    mutate(RankMinSD=length(MinSD)-rank(MinSD, ties.method="first")) %>%
    ungroup()

pl <- sd_min %>%
    dplyr::filter(Group=="lnoncoding") %>% 
    ggplot(., aes(x=MinSD, y=RankMinSD)) +
    geom_line(size=1) +
    labs(y="#Genes", x="Minimum Standard Deviation", title="Gene Filtering\nLong Noncoding RNA") +
    xlim(0, 1.5) +
    theme_minimal()
ggsave(pl, filename=file.path("figures", "filter_sd_cemitool", "lnoncoding_minsd.pdf"))



