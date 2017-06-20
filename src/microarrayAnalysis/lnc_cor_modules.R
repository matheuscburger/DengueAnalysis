#!/usr/bin/env Rscript

"Calculate correlation between lncRNAs and genes in modules from CEMiTool

Usage: lnc_cor_modules.R --modules=<file> --expression=<file> --gencode=<file> --biotype=<file> --output-corr=<file> --output-pvalue=<file> --output-padj=<file> --output-stats=<file> [--cutoff-cor=<value> --cutoff-pvalue=<value> --cutoff-padj=<value> --ensembl-col=<value> --annotation-cols=<value>...]

Options:
  -h --help                   show this help message
  --version                   show program version
  --modules=<file>            a file containing the gene modules (GMT format)
  --expression=<file>         a expression file (format: TSV)
  --gencode=<file>            gencode annotation, should have ensembl gene ids and biotype name
  --biotype=<file>            a file from Gencode database, should have object_type, name and biotype_group
  --output-corr=<file>        output file to write correlations between lncRNAs and genes in all modules
  --output-pvalue=<file>      output file to write p-values calculated from correlation values
  --output-padj=<file>        output file to write adjusted p-values
  --output-stats=<file>       output file to write results
  --cutoff-cor=<value>        correlation cutoff [default: 0.7]
  --cutoff-pvalue=<value>     p-value cutoff [default: 0.05]
  --cutoff-padj=<value>       adjusted p-value [default: 0.1]
  --ensembl-col=<value>       column containing ENSEMBL ID [default: Symbol]
  --annotation-cols=<value>   annotation columns

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

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


if (!interactive()) {
    # Get and check arguments.
    suppressMessages(library("docopt"))

    arg <- docopt(doc, version="0.0.1\n", strict=T)
    arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
    clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
    names(arg) <- clean(names(arg))

    parameters <- arg
    modules_fname <- parameters[["modules"]]
    exp_fname <- parameters[["expression"]]
    gencode_annot_fname <- parameters[["gencode"]]
    biotypes_fname <- parameters[["biotype"]]
    annotation_cols <- parameters[["annotation_cols"]]
    ensembl_col <- parameters[["ensembl_col"]]
    output_cor_fname <- parameters[["output_corr"]]
    output_ps_fname <- parameters[["output_pvalue"]]
    output_padj_fname <- parameters[["output_padj"]]
    output_stats_fname <- parameters[["output_stats"]]
    cutoff_cor <- as.numeric(parameters[["cutoff_cor"]])
    cutoff_p <- as.numeric(parameters[["cutoff_pvalue"]])
    cutoff_p_adj <- as.numeric(parameters[["cutoff_padj"]])
    annotation_cols <- c(annotation_cols, ensembl_col)

    library("readr")
    library("dplyr")
    library("tidyr")
    library("stringr")

    # get gene annotation
    message("Processing gene annotations ...")
    biotypes <- read_tsv(biotypes_fname) %>% 
        dplyr::filter(object_type == "gene") %>% 
        dplyr::select(name, biotype_group)
    gencode_annot <- read_tsv(gencode_annot_fname, col_names=c("ENS_ID", "name")) %>%
        mutate(ENS_ID=str_replace(ENS_ID, "\\.\\d+$", "")) %>%
        left_join(biotypes, by=c("name"="name"))

    group2gene <- split(gencode_annot[["ENS_ID"]], gencode_annot[["biotype_group"]])

    # get gene expression
    message("Reading gene expression table ...")
    exp_tib <- read_tsv(exp_fname)
    exp_mat <- exp_tib %>% select_(paste0("-", annotation_cols)) %>% as.matrix(.)
    rownames(exp_mat) <- exp_tib[[ensembl_col]]

    # get modules
    message("Reading modules (gmt) file ...")
    modules <- read.gmt(modules_fname)[["genes"]]

    all_modules_genes <- unlist(modules)
    names(all_modules_genes) <- NULL


    # extract module genes from gene expression matrix
    message("Subseting genes expression matrix ...")
    genes_to_extract <- rownames(exp_mat)[which(rownames(exp_mat) %in% all_modules_genes)]
    mod_mat <- exp_mat[genes_to_extract, ]
    rm(genes_to_extract)

    # calculate mean expression of each module
    message("Calculating mean expression for each module ...")
    module_mean <- lapply(modules, function(x) colMeans(mod_mat[x, ], na.rm=T)) %>%
        do.call(rbind, .)

    # extract lnoncoding from gene expression matrix
    message("Subseting lncRNAs from gene expression matrix ...")
    genes_to_extract <- rownames(exp_mat)[which(rownames(exp_mat) %in% group2gene[["lnoncoding"]])]
    lnc_mat <- exp_mat[genes_to_extract, ]
    rm(genes_to_extract)

    # calculate correlation
    message("Transposing matrices ...")
    lnc_transposed <- t(lnc_mat)
    mod_transposed <- t(mod_mat)
    mod_mean_transposed <- t(module_mean)

    message("Calculating correlation between genes in modules and lncRNAs ...")
    cor_lnc_mod <- cor(lnc_transposed, mod_transposed)
    write.table(cor_lnc_mod, output_cor_fname, sep="\t", col.names=NA, row.names=T, quote=FALSE)

    message("Calculating p-values ...")
    ps <- calc.p(cor_lnc_mod, ncol(exp_mat))
    write.table(ps, output_ps_fname, sep="\t", col.names=NA, row.names=T, quote=FALSE)

    message("Calculating adjusted p-values ...")
    adj_ps <- matrix(p.adjust(as.vector(ps), method="fdr"), ncol=ncol(ps))
    write.table(ps, output_padj_fname, sep="\t", col.names=NA, row.names=T, quote=FALSE)

    message("Calculating correlation between mean module expression and lncRNAs ...")
    cor_lnc_mean <- cor(lnc_transposed, mod_mean_transposed)

    # 
    list_res <- list()

    message("Subseting correlation matrices ...")
    cor_bool <- abs(cor_lnc_mod) > cutoff_cor
    p_bool <- ps < cutoff_p
    padj_bool <- adj_ps < cutoff_p_adj
    cor_above_cutoff <- cor_bool & p_bool & padj_bool
    for(mod in names(modules)){
        message("Processing ", mod, " ...")
        mod_genes <- modules[[mod]]
        cor_above_cutoff_mod <- cor_above_cutoff[,mod_genes]

        # mean correlation with all genes in module
        mean_cor <- apply(cor_lnc_mod[, mod_genes, drop=F], 1, mean)
        
        # count the number of genes in the module correlated with a lncRNA
        count <- sort(rowSums(cor_above_cutoff_mod), decreasing=T)
        count_0 <- which(count == 0)
        if(length(count_0) > 0){
            count <- count[-count_0]
        }

        if(length(count) > 0){
            message(mod, " has at least one lncRNA correlated.")
            # which genes
            genes <- split(cor_above_cutoff_mod[names(count), ], 1:length(count)) %>%
                lapply(., function(x) mod_genes[which(x)]) %>%
                lapply(., paste, collapse=",") %>%
                unlist

            l <- data_frame(lncRNA=names(count), module=mod, N_Correlated=count, 
                            N_Module=length(mod_genes), Proportion=count/length(mod_genes),
                            corr_with_mean=cor_lnc_mean[names(count), mod],
                            mean_correlation_all_genes=mean_cor[names(count)],
                            genes=genes)
            list_res[[mod]] <- l
        }
    }

    message("Joining results for each module... ")
    res <- do.call(rbind, list_res)

    message("Writing output ...")
    write_tsv(res, output_stats_fname)
}
