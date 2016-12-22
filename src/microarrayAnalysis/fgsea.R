#!/usr/bin/env Rscript

"FGSEA

Usage: fgsea.R --input=<file> --output=<file> --gmt=<file> --symbols=<value> [--annotation-cols=<value>...]

Options:
  -h --help                  show this help message
  --version                  show program version
  --input=<file>             input file
  --output=<file>            output file
  --gmt=<file>               gmt file
  --symbols=<value>          column containing symbols
  --class-col=<value>        column containing class in sample annotation file
  --annotation-cols=<value>  column containing gene annotations

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

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
    res[['genes']] <- gmt.genes
    res[['desc']] <- gmt.desc
    return(res)
}

onesample_fgsea <- function(ranks, gene_symbols){
	names(ranks) <- gene_symbols
	fgsea_res <- fgsea(pathways = geneset, 
					  stats = ranks,
					  minSize=15,
					  maxSize=500,
					  nperm=10000,
					  nproc=0,
					  )
	fgsea_res[['leadingEdge']] <- sapply(fgsea_res[['leadingEdge']], paste, collapse=",")
	return(fgsea_res)
}

run_fgsea <- function(input_gsea, genesets, symbols_col){

	register(SerialParam())
	bpparameters <- bpparam() 

	gene_symbols <- input_gsea[[symbols_col]]

	input_gsea <- input_gsea %>% select_(paste0("-",symbols_col))

	scores <- input_gsea %>% lapply(., onesample_fgsea, gene_symbols)
	for(n in names(scores)){
	    colnames(scores[[n]])[-1] <- paste0(n, "_", colnames(scores[[n]])[-1])
	}
	allscores <- Reduce(function(x, y) merge(x, y, by="pathway", all=TRUE), scores)
        
    return(allscores)
}    

suppressMessages(library('fgsea'))
suppressMessages(library('dplyr'))
suppressMessages(library('readr'))
suppressMessages(library('data.table'))
suppressMessages(library('BiocParallel'))

if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	# parameters

	message("Reading TSV ...")
	input.tibble <- read_tsv(arg$input) %>%
		filter_(paste0("!is.na(", arg$symbols,")"))

	if("annotation_cols" %in% names(arg)){
		message("Removing annotation columns ...")
		input.tibble <- input.tibble %>% 
			dplyr::select_(.dots=paste0("-", arg$annotation_cols))
	}

	message("Reading GMT files ...")
	gmt <- read.gmt(arg$gmt)
	geneset <- gmt[["genes"]]

	message("Running run_fgsea ...")
	res <- run_fgsea(input.tibble, geneset, arg$symbols)

	message("Writing output file ...")
	write_tsv(x=res, path=arg$output)
}
