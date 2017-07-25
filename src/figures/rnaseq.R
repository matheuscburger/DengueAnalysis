#!/usr/bin/env Rscript

library("tidyverse")
library("forcats")
library("data.table")
library("RColorBrewer")
library("ggrepel")
library("gridExtra")
library("gtable")
library("grid")
library("stringr")
library("scales")

dir.create("figures/rnaseq")

# function to fix the reannotation names in order to match rnaseq names
fix_reannot_names <- function(gene, database){
    .fix_reannot_names <- function(gene_id, database){
        #gene_id <- x['Gene']
        #database <- x['Database']
        if (database == "gencode"){
            # do nothing
        } else if (database == "noncode"){
            # do nothing
        } else if (database == "lncipedia") {
            gene_id <- sub(";_alignment_id_\\d+$", "", gene_id)
        } else if (database == "mitranscriptome"){
            gene_id <- sub("^gene_id_", "", gene_id)
        }
        return(gene_id)
    }
    new_gene <- character(length(gene))
    for(i in 1:length(gene)){
        new_gene[i] <- .fix_reannot_names(gene[i], database[i])
    }
    return(new_gene)
}
# ler tabela de anotacao
plats <- c("GPL570", "GPL2700")
annotation_fname <- "config/reannotation/annotation_long.tsv"
reannot <- read_tsv(annotation_fname) %>% 
    filter(Hits == 1, NumAnnot == 1, Platform %in% plats) %>%
    dplyr::select(Gene, Platform, Database, Probe, Type, Group) %>%
    mutate(Gene=fix_reannot_names(Gene, Database)) %>%
    group_by(Gene, Platform) %>%
    mutate(Probe=paste(sort(Probe), collapse="/")) %>%
    unique()  %>%
    spread(key=Platform, value=Probe) %>%
    ungroup()

# ler dados de rnaseq
rnaseq_dir <- "config/rnaseq_data/genes_pme_TPM"

rnaseq_files <- grep("GSE52166", list.files(rnaseq_dir, full.names=T), value=T)


rnaseq_dt <- lapply(rnaseq_files, fread)
names(rnaseq_dt) <- gsub("\\.tsv", "", basename(rnaseq_files))
# include the name of study in sample columns
for( n in names(rnaseq_dt) ){
    names(rnaseq_dt[[n]])[-1] <- paste(n, names(rnaseq_dt[[n]]), sep="_")[-1]
}

# ler tabela de outliers
outliers <- readLines("config/rnaseq_data/outliers.txt")
out <- paste(outliers, collapse="|")

# incluir gene types
gencode2type <- read_tsv("config/reannotation/gencode_annotation.tsv", col_names=c("Gene", "Type")) %>%
    mutate(Gene=str_match(Gene, "^([A-Z]+\\d+)\\.\\d+")[,2]) %>%
    with(., split(Type, Gene)) %>%
    as.environment()

# pegar groups de cada gene types
type2group <- read_tsv("config/reannotation/biotypes.tsv") %>%
    filter(is_current==1) %>%
    select(name, biotype_group) %>%
    unique() %>%
    with(., split(biotype_group, name)) %>%
    as.environment()


get_database <- function(x){
    sapply(x, function(y){
        if (grepl("^ENSG\\d+", y)){
            return("gencode")
        } else if (grepl("NONHSAG\\d+", y)){
            return("noncode")
        } else if (grepl("^lnc-", y)){
            return("lncipedia")
        } else if (grepl("G\\d+", y)){
            return("mitranscriptome")
        } else {
            return(NA)
        }
    })
}

# juntar estudos e remover outliers
rnaseq_merged <- as_data_frame(reduce(rnaseq_dt, merge, by="gene_id")) %>%  
    mutate(Gene=str_match(gene_id, "^([A-Z]+\\d+)\\.\\d+")[,2]) %>%
    mutate(gene_id=ifelse(is.na(Gene), gene_id, Gene)) %>%
    select(-Gene, -matches(out)) %>%
    rowwise() %>% 
    mutate(Database=get_database(gene_id),           
           Type=unlist(mget(gene_id, gencode2type, ifnotfound=NA))) %>% 
    mutate(Group=unlist(mget(Type, type2group, ifnotfound=NA))) %>%
    mutate(Type=ifelse(Database=="gencode", Type, "noncoding"),
           Group=ifelse(Database=="gencode", Group, "lnoncoding")) %>% 
    right_join(select(reannot, Gene, starts_with("GPL")), .,  by=c("Gene"="gene_id")) %>%
    select(-starts_with("GSE"), starts_with("GSE"))


rank_decreasing <- function(x) rank(-x)

rnaseq_rank <- rnaseq_merged %>% 
    ungroup() %>%
    mutate_at(vars(starts_with("GSE")), funs(rank_decreasing))


# figuras

# figura 1

# juntar os estudos ? ok 
# calcular Rank - ok 
# pegar genes que aparecem em pelo menos 90% das amostras
# entre os 1000, 2000, 3000, 4000, ... genes com maior expressão

rnaseq_melt <- rnaseq_rank %>% 
    filter(!(is.na(GPL570) & is.na(GPL2700))) %>% 
    gather("Study_Sample", "Expression", starts_with("GSE")) %>%
    separate(Study_Sample, c("Study", "Sample"), sep="_")

#topgenes <- rnaseq_melt %>%  mutate(Top1000=Expression < 1e3,
#                                    Top2000=Expression < 2e3,
#                                    Top3000=Expression < 3e3,
#                                    Top4000=Expression < 4e3)

gen_top <- function(x, cuts){
    top_list <- lapply(cuts, function(y) x < y)
    names(top_list) <- paste0("Top_", cuts)
    return(as_data_frame(top_list))
}
cuts <- seq(1e3, 2e4, by=1e3)
topgenes <- rnaseq_melt %>% ungroup() %>% 
    do(gen_top(.[["Expression"]], cuts)) %>%
    bind_cols(rnaseq_melt, .)

nsamples <- length(unique(topgenes[["Sample"]]))

group_levels <- c("coding", "unannotated", "lnoncoding",
                  "pseudogene", "snoncoding", "mnoncoding",
                  "undefined")
change_colors <- function(incolors){
    myhsv <- rgb2hsv(col2rgb(incolors))
    myhsv["s", ] <- sapply(myhsv["s", ] + 0.8, min, 1)
    myhsv["v", ] <- sapply(myhsv["v", ] - 0.4, max, 0)
    res <- apply(myhsv, 2, function(x) hsv(h=x[1], s=x[2], v=x[3]))
    return(res)
}
group_colors <- change_colors(brewer.pal(length(group_levels), "Set1"))
names(group_colors) <- group_levels

sample_freq <- function(x) {(sum(x, na.rm=T)/nsamples)>.9}
top_counts <- topgenes %>% select(Gene, Group, starts_with("Top")) %>% 
    ungroup %>% group_by(Gene) %>% 
    mutate_at(vars(starts_with("Top")), funs(sample_freq)) %>%
    unique() %>% 
    ungroup() %>% group_by(Group) %>%
    summarise_at(vars(starts_with("Top")), funs(sum)) %>%
    gather("Top", "Count", -Group) %>% 
    mutate(Top=as_factor(Top, paste0("Top_", cuts)),
           Group=as_factor(Group, levels=group_levels))

top_sum <- top_counts %>% group_by(Top) %>% summarise(Sum=sum(Count))

# fazer um gráfico de barras mostrando os types e group_types 
# simplifica tema para grafico
p1 <- ggplot(top_sum, aes(x=Top, y=Sum)) +
    geom_bar(stat="identity", fill="#E59C00") +
    geom_text(aes(label=Sum), vjust=0, nudge_y=50, color="#304D78", fontface="bold", size=5) +
    labs(title="Microarray Genes Among Most Expressed Genes") + 
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())

p2 <- ggplot(filter(top_counts, Group=="lnoncoding"), aes(x=Top, y=Count, group=1)) +
    geom_line() +
    geom_point() +
    geom_label(aes(label=Count), vjust=0, nudge_y=5, fontface="bold", size=5, label.size=0, alpha=.8, fill="white") +
    labs(title="Number of Long non-coding RNAs", y="LncRNA") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1))


g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
# avoid geom_text clipping
g1$layout$clip[which(g1$layout$name == "panel")] <- "off"
g2$layout$clip[which(g2$layout$name == "panel")] <- "off"
colnames(g1) <- paste0(seq_len(ncol(g1)))
colnames(g2) <- paste0(seq_len(ncol(g2)))
pdf("figures/rnaseq/barplot_topX.pdf", width=10)
grid.newpage()
grid.draw(combine(g1, g2, along=2))
dev.off()

# figura 2 

# gráfico de densidade para cada type e group_type
rnaseq_mean <- rnaseq_merged %>% 
    mutate(Mean=rowMeans(select(., starts_with("GSE52166")))) %>%
    select(-starts_with("GSE52166"))

tpm_cut <- 10
pl <- rnaseq_mean %>%
    filter(Group %in% c("coding", "lnoncoding")) %>%
    ggplot(., aes(x=Mean, fill=Group)) + 
        geom_density(alpha=.5) +
        scale_fill_manual(values=group_colors) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        geom_vline(xintercept=tpm_cut) +
        annotate("text", x=tpm_cut, y=Inf, label=paste0("x = ", tpm_cut), vjust=2, hjust=-0.1, fontface="bold") +
        annotation_logticks(sides="trbl") +
        labs(title="Density of Mean Expression (Log10 Scale)\nAll Genes", x="Mean TPM") +
        theme_bw()
ggsave(plot=pl, filename="figures/rnaseq/density_mean_TPM.pdf")


pl <- rnaseq_mean %>%
    filter(Group %in% c("coding", "lnoncoding"),
           !(is.na(GPL570) & is.na(GPL2700))) %>%
    ggplot(., aes(x=Mean, fill=Group)) + 
        geom_density(alpha=.5) +
        scale_fill_manual(values=group_colors) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        geom_vline(xintercept=tpm_cut) +
        annotate("text", x=tpm_cut, y=Inf, label=paste0("x = ", tpm_cut), vjust=2, hjust=-0.1, fontface="bold") +
        annotation_logticks(sides="trbl") +
        labs(title="Density of Mean Expression (Log10 Scale)\nGenes Containing Reliable Probes", x="Mean TPM") +
        theme_bw()
ggsave(plot=pl, filename="figures/rnaseq/density_mean_TPM_probes.pdf")

# figura 3
# definir um corte de TPM com base no gráfico de densidade
# fazer um gráfico de barras com os types e group_types
tpm_count_after_cut <- rnaseq_mean %>% filter(Mean >= tpm_cut) %>%
    count(Group)

pl <- ggplot(tpm_count_after_cut, aes(x=Group, y=n, fill=Group)) +
    geom_bar(stat="identity") +
    geom_text(aes(y=n/2, label=n), color="white", fontface="bold") +
    scale_fill_manual(values=group_colors) +
    labs(title=paste0("Number of Genes (TPM >= ", tpm_cut, ")\n",
                      "All Genes"), y="Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(filename="figures/rnaseq/barplot_after_cut.pdf", plot=pl)

tpm_count_after_cut_probes <- rnaseq_mean %>% 
    filter(Mean >= tpm_cut, 
           !(is.na(GPL570) & is.na(GPL2700))) %>%
    count(Group)

pl <- ggplot(tpm_count_after_cut_probes, aes(x=Group, y=n, fill=Group)) +
    geom_bar(stat="identity") +
    geom_text(aes(y=n/2, label=n), color="white", fontface="bold") +
    scale_fill_manual(values=group_colors) +
    labs(title=paste0("Number of Genes (TPM >= ", tpm_cut, ")\n",
                      "Genes Containing Reliable Probes"), y="Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(filename="figures/rnaseq/barplot_after_cut_probes.pdf", plot=pl)
