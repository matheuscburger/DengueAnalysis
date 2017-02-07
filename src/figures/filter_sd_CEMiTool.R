library("tidyverse")
library("stringr")

dir.create(file.path('figures', 'filter_sd_cemitool'))

expression_files_directory <- "/home/mburger/work/mburger/filter_cemitool/collapsed"
reannotation_filename <- "/home/mburger/work/mburger/filter_cemitool/reannotation/annotation_long.tsv"

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
