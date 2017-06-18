

library("readr")
library("dplyr")
library("tidyr")
library("stringr")


calc.p <- function(r, n){ 
    r <- abs(r) 
    r.1 <- matrix(mapply(function(x, y) {isTRUE(all.equal(x, y))}, r, 1),  
           nrow=nrow(r), ncol=ncol(r)) 
    r[which(r.1, arr.ind=T)] <- 1 
    t.val <- ( r * sqrt(n-2) ) / sqrt(1 - r**2) 
    p <- pt(t.val, df=n-2, lower.tail=FALSE)*2 
    return(p) 
}

p2r <- function(p, n){
    dfreedom <- n-2
    t.val <- qt(p/2, df=dfreedom, lower.tail=F)
    return(t.val/sqrt(t.val**2+dfreedom))
}

read.gmt <- function(fname){
    res <- list(genes=list(), desc=list())
    gmt <- file(fname)
    gmt.lines <- readLines(gmt)
    close(gmt)
    gmt.list <- lapply(gmt.lines, function(x) unlist(strsplit(x, split="\t")))
    gmt.names <- sapply(gmt.list, '[', 1)
    gmt.desc <- lapply(gmt.list, '[', 2)
    gmt.genes <- lapply(gmt.list, function(x){x[3:length(x)]})
    names(gmt.desc) <- names(gmt.genes) <- gmt.names
    return(list(genes=gmt.genes, desc=gmt.desc))
}

modules_fname <- "results/CEMiTool_joined_modules.gmt"
exp_fname <- "data/processed/collapsed/GSE43777.tsv"
gencode_annot_fname <- "config/reannotation/gencode_annotation.tsv"
biotypes_fname <- "config/reannotation/biotypes.tsv"
annotation_cols <- c()
ensembl_col <- "Symbol"
output_cor_fname <- "results/CEMiTool_joined/corr_lnc/GSE43777.tsv"
output_ps_fname <- "results/CEMiTool_joined/corr_pvalue_lnc/GSE43777.tsv"
output_padj_fname <- "results/CEMiTool_joined/corr_padj_lnc/GSE43777.tsv"
output_stats_fname <- "results/CEMiTool_joined/corr_stats_lnc/GSE43777.tsv"
cutoff_cor <- 0.8
cutoff_p <- 1
cutoff_p_adj <- 1

annotation_cols <- c(annotation_cols, ensembl_col)

# get gene annotation
biotypes <- read_tsv(biotypes_fname) %>% 
    dplyr::filter(object_type == "gene") %>% 
    dplyr::select(name, biotype_group)
gencode_annot <- read_tsv(gencode_annot_fname, col_names=c("ENS_ID", "name")) %>%
    mutate(ENS_ID=str_replace(ENS_ID, "\\.\\d+$", "")) %>%
    left_join(biotypes, by=c("name"="name"))

group2gene <- split(gencode_annot[["ENS_ID"]], gencode_annot[["biotype_group"]])

# get gene expression
exp_tib <- read_tsv(exp_fname)
exp_mat <- exp_tib %>% select_(paste0("-", annotation_cols)) %>% as.matrix(.)
rownames(exp_mat) <- exp_tib[[ensembl_col]]

# get modules
modules <- read.gmt(modules_fname)[["genes"]]

all_modules_genes <- unlist(modules)
names(all_modules_genes) <- NULL


# extract module genes from gene expression matrix
genes_to_extract <- rownames(exp_mat)[which(rownames(exp_mat) %in% all_modules_genes)]
mod_mat <- exp_mat[genes_to_extract, ]
rm(genes_to_extract)

# extract lnoncoding from gene expression matrix
genes_to_extract <- rownames(exp_mat)[which(rownames(exp_mat) %in% group2gene[["lnoncoding"]])]
lnc_mat <- exp_mat[genes_to_extract, ]
rm(genes_to_extract)

# calculate correlation
lnc_transposed <- t(lnc_mat)
mod_transposed <- t(mod_mat)

cor_lnc_mod <- cor(lnc_transposed, mod_transposed)
write.table(cor_lnc_mod, output_cor_fname, sep="\t", col.names=NA, row.names=T, quote=FALSE)

ps <- calc.p(cor_lnc_mod, ncol(exp_mat))
write.table(ps, output_ps_fname, sep="\t", col.names=NA, row.names=T, quote=FALSE)

adj_ps <- matrix(p.adjust(as.vector(ps), method="fdr"), ncol=ncol(ps))
write.table(ps, output_padj_fname, sep="\t", col.names=NA, row.names=T, quote=FALSE)

# 
list_res <- list()

cor_bool <- abs(cor_lnc_mod) > cutoff_cor
p_bool <- ps > cutoff_p
padj_bool <- adj_ps > cutoff_p_adj
cor_above_cutoff <- cor_bool & p_bool & padj_bool
for(mod in names(modules)){
    mod_genes <- modules[[mod]]
    cor_above_cutoff_mod <- cor_above_cutoff[,mod_genes]
    
    # count the number of genes in the module correlated with a lncRNA
    count <- sort(rowSums(cor_above_cutoff_mod), decreasing=T)
    count_0 <- which(count == 0)
    if(length(count_0) > 0){
        count <- count[-count_0]
    }

    if(length(count) > 0){
        # which genes
        genes <- split(cor_above_cutoff_mod[names(count), ], 1:length(count)) %>%
            lapply(., function(x) mod_genes[which(x)]) %>%
            lapply(., paste, collapse=",") %>%
            unlist

        l <- data_frame(lncRNA=names(count), module=mod, n=count, genes=genes)
        list_res[[mod]] <- l
    }
}

res <- do.call(rbind, list_res)

write_tsv(res, output_stats_fname)
