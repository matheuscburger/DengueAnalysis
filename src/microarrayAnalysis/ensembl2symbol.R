#!/usr/bin/env Rscript

"Ensembl ID To Gene Symbols

Usage: ensembl2symbol.R --input=<file> --output=<file> [--ensembl-col=<value> --is-list --attributes=<value>...] 

Options:
  -h --help                  show this help message
  --version                  show program version
  --input=<file>             input file name (table or list)
  --output=<file>            output file name
  --ensembl-col=<value>      column containing ensembl ids
  --attributes=<value>       aditional attributes [default: external_gene_name]
  --is-list                   input file is a list

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc


suppressMessages(library("readr"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("stringr"))
suppressMessages(library("biomaRt"))

if (!interactive() && !exists('SOURCE')) {
    # Get and check arguments.
    suppressMessages(library(docopt))
    arg <- docopt(doc, version="0.0.1\n", strict=T)
    arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
    clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
    names(arg) <- clean(names(arg))


    ## RUN
    #arg <- list()
    #arg$input <- "data/processed/filtered/GSE13052.tsv"
    #arg$ensembl_col <- "Symbol"
    #arg$attributes <- c("external_gene_name")
    #arg$output <- "comsymbols.tsv"

    if("is_list" %in% names(arg)){
        ensgs <- sort(unique(readLines(arg$input)))
    } else {
        in_df <- read_tsv(arg$input)
        ensgs <- sort(unique(in_df[[arg$ensembl_col]]))
    }

    ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="ensembl.org")
    gsymbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", arg$attributes), filters="ensembl_gene_id",
                      values=ensgs, mart=ensembl)


    if(!arg[["is_list"]]){
        by_vec <- arg$ensembl_col
        names(by_vec) <- "ensembl_gene_id"
        gsymbols <- in_df %>% right_join(gsymbols, ., by=by_vec) 
    }
    gsymbols %>% write_tsv(arg$output)
}
