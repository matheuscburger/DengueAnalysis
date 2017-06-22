#!/usr/bin/env Rscript

"join_lnc_cor_modules.R

Usage: join_lnc_cor_modules.R --output=<file> [--join-by=<value>... --genes-suffix=<value> --n-correlated-col=<value>] INPUT...

Options:
  -h --help                    show this help message
  --version                    show program version
  INPUT                        input file names
  --output=<file>              output file name
  --n-correlated-col=<value>   Column containing the number of genes correlated. [default: N_Correlated]
  --join-by=<value>            Columns to use in merge process. [default: lncRNA, module]
  --genes-suffix=<value>       Suffix used in columns containing genes. [default: _genes]

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

if (!interactive() && !exists('SOURCE')) {
    # Get and check arguments.
    suppressMessages(library("docopt"))
    arg <- docopt(doc, version="0.0.1\n", strict=T)
    arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
    clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
    names(arg) <- clean(names(arg))

    library("data.table")
    library("stringr")
    library("tibble")
    library("dplyr")
    library("readr")

    input_fnames <- arg[["input"]]
    output_fname <- arg[["output"]]
    join_by <- arg[["join_by"]]
    if(join_by == "lncRNA, module"){
        join_by <- c("lncRNA", "module")
    }
    genes_suffix <- arg[["genes_suffix"]]
    n_correlated_col <- arg[["n_correlated_col"]]

    number_inputs <- length(input_fnames)

    files <- lapply(input_fnames, fread)
    names(files) <- str_replace(basename(input_fnames), "\\.tsv$", "")

    for(n in names(files)){
        cols_idx <- which(!colnames(files[[n]]) %in% join_by)
        colnames(files[[n]])[cols_idx] <- paste0(n, "_", colnames(files[[n]])[cols_idx])
    }

    mrgd <- as_data_frame(Reduce(function(x, y) merge(x, y, by=join_by, all=T), files)) %>%
        mutate(N_Studies=rowSums(!is.na(dplyr::select(., ends_with(n_correlated_col))))) %>%
        mutate(Prop_Studies=N_Studies/number_inputs)

    genes_cols <- paste0(names(files), genes_suffix)
    intersect_genes <- apply(mrgd[, genes_cols], 1, function(x) Reduce(intersect, strsplit(x[which(!is.na(x))], ",")))
    n_intersect <- sapply(intersect_genes, length)
    paste_int_genes <- sapply(intersect_genes, paste0, collapse=",")

    mrgd <- mrgd %>% 
        mutate(intersect_genes=paste_int_genes, N_intersect=n_intersect) %>%
        arrange(desc(N_Studies), desc(N_intersect))

    write_tsv(mrgd, output_fname)
}
