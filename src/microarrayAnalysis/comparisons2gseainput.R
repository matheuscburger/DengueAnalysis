#!/usr/bin/env Rscript

"Get input for GSEA from statistics file.

Usage: comparisons2gseainput.R --input=<file> --output=<file> [--ensembl-col=<value> --stats-col=<value> --comparison-col=<value> --dontconvert]

Options:
  -h --help                  show this help message
  --version                  show program version
  --input=<file>             input file name
  --output=<file>            output file name
  --ensembl-col=<value>      column containing ensembl ids [default: Symbol]
  --stats-col=<value>        column containing statistics [default: Log2FC]
  --comparison-col=<value>   column containing comparison title [default: title]
  --dontconvert              don't convert ensembl ids to symbols

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

library("ggplot2")
library("dplyr")
library("tidyr")
library("tibble")
library("readr")
library("biomaRt")

if (!interactive() && !exists('SOURCE')) {
    # Get and check arguments.
    suppressMessages(library(docopt))
    arg <- docopt(doc, version="0.0.1\n", strict=T)
    arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
    clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
    names(arg) <- clean(names(arg))


    ## RUN
    fname <- arg$input
    stats.df <- read_tsv(fname)

    stats_column <- arg$stats_col
    gene_column <- arg$ensembl_col
    comparison_column <- arg$comparison_col

    wide <- stats.df %>% 
        dplyr::select_(gene_column, comparison_column, stats_column) %>% # select columns of interest
        .[!is.na(.[[gene_column]]), ] %>%  # remove NAs
        group_by_(gene_column, comparison_column) %>%
        filter_(paste0("abs(", stats_column, ") == max(abs(", stats_column,"))")) %>% 
        slice(1) %>% ungroup() %>% 
        spread_(key=comparison_column, value=stats_column)


    if (!arg[["dontconvert"]]){
        message("Trying to convert Ensembl IDs to Symbols...")
        ensgs <- wide[[gene_column]]

        ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="ensembl.org")
        gsymbols <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters="ensembl_gene_id",
                          values=ensgs, mart=ensembl)

        wide <- wide %>% 
            right_join(as_tibble(gsymbols), ., by=c('ensembl_gene_id'=gene_column)) %>% 
            dplyr::select(-ensembl_gene_id)
    }
    write_tsv(wide, arg$output)
}
