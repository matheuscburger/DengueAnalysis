#!/usr/bin/env Rscript

library("readr")
library("stringr") 
library("dplyr")
library("tidyr")
library("corrplot")
library("forcats")
library("ggplot2")
library("scales")

dir.create("figures/CEMiTool_joined_corrplot")

fgsea_files <- list.files("results/CEMiTool_joined/fgsea", full.names=T)

fgsea_l <- lapply(fgsea_files, read_tsv)

names(fgsea_l) <- str_extract(basename(fgsea_files), "GSE\\d+")

for(n in names(fgsea_l)){
    fgsea_l[[n]] <- fgsea_l[[n]] %>% mutate(Study=n)
}

nes <- fgsea_l %>% lapply(., function(x) x %>%
                          select(Study, pathway, ends_with("_NES")) %>%
                          gather(Comparison, NES, -Study, -pathway)) %>%
    do.call(rbind, .) %>%
    mutate(Comparison=sub("_NES", "", Comparison))


pvalue <- fgsea_l %>% lapply(., function(x) x %>%
                          select(Study, pathway, ends_with("_padj")) %>%
                          gather(Comparison, Adj_PValue, -Study, -pathway)) %>%
    do.call(rbind, .) %>%
    mutate(Comparison=sub("_padj", "", Comparison))

fgsea_tidy <- inner_join(nes, pvalue, by=c("Study", "pathway", "Comparison"))

fgsea_filter <- fgsea_tidy %>%
    filter(Adj_PValue <= 0.05)

my_squish <- function(...){
    return(scales::squish(..., only.finite=FALSE))
}

corrColors <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
pl <- fgsea_filter %>% 
    group_by(Comparison, pathway) %>%
    mutate(Sum=sum(NES)) %>%
    group_by(Comparison) %>%
    arrange(Sum) %>%
    do(pl = ggplot(., aes(y=fct_inorder(pathway), x=Study, fill=NES, size=abs(NES))) +
            geom_point(color="white", shape=21) +
            scale_fill_gradientn(colours=corrColors,
                                 limits=c(-2, 2),
                                 oob=my_squish) +
            scale_size_area(limits=c(0, 2), max_size=20, oob=my_squish, guide="none") +
            scale_x_discrete(position = "top") +
            labs(title=unique(.[["Comparison"]]),
                 x="Study", y="Module") +
            theme_classic(base_size=18) +
            theme(panel.background = element_blank(),
                  axis.line = element_blank())
        )


pdf("figures/CEMiTool_joined_corrplot/NES.pdf", height=12, width=7)
.null <- lapply(pl["pl"], print)
dev.off()
