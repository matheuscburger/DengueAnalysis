#!/usr/bin/env Rscript

"Annotate probes in a table

Usage: annotate_probes.R INPUT OUTPUT (--annotation-package=<pack> | --annotation-file=<file>) [--multiple=<value>]

Input:
  INPUT                      File to annotate (Tab Separated Values format), first column must have probe ID

Output:
  OUTPUT                     Annotated file

Options:
  -h --help                    show this help message
  --version                    show program version
  --annotation-file=<file>     annotation file containing the probe id in the first column e symbol gene in the second column
  --annotation-package=<pack>  annotation package (AnnotationDbi)
  --multiple=<value>           What should I do with multiple genes mapping to one probe? [default: join]


Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

#' Get a annotation environment that converts probe in gene symbol
#'
#' @param fname	file with mapping
#' @return	an environment
get_annot_env <- function(fname){
	probe_annot <- read.delim(fname, stringsAsFactors=FALSE, header=FALSE)
	tmp_list <- as.list(tapply(probe_annot[,2], probe_annot[,1], c))
	return(as.environment(tmp_list))
}


#' Add a column with Gene Symbols in a data.frame.
#'
#' @param exp_df	a data.frame with expression values for each probe (rownames are probenames)
#' @param annot_env	a environment that converts probe name in gene symbol
#' @param probe_col column containing probe names
#' @param remove_pms some affy platforms have a suffix "_PM_", should I remove this before the annotation?
#' @param multiple  what should I do with multiple symbols for one probe (NA, first, first_sorted, join)
#' @return a data frame with one new column called Symbol
annotate <- function(exp_df, annot_env, probe_col="row.names", remove_pms=TRUE, multiple=c("join", "first", "first_sorted", "NA"), ...){
	if(multiple == "NA"){
		multiple <- "mNA"
	}
	join <- function(x){
		return(paste(x, sep="/", collapse="/"))
	}
	mNA <- function(x){
		if(length(x) > 1){
			return(NA)
		}else{
			return(x)
		}
	}
	first <- function(x){
		return(x[1])
	}
	first_sorted <- function(x){
		return(sort(x)[1])
	}

	if(probe_col == "row.names"){
		pnames <- as.character(rownames(exp_df))
	}else{
		pnames <- as.character(exp_df[, probe_col])
	}

	if(remove_pms && any(grep("_PM_", pnames))){
		pnames <- sub("_PM_", "_", pnames)
	}

	symbols <- mget(pnames, annot_env, ifnotfound=NA)
	ssymbols <- sapply(symbols, get(multiple))
	exp_df <- cbind(Symbol=ssymbols, exp_df)

	return(exp_df)
}

annotate_main <- function(input, output, annotation_package, annot_env, probe_col, multiple, ...){
	if(!missing(annotation_package)){
		suppressMessages(library(paste(annotation_package), character.only=TRUE))
		annot_env <- toggleProbes(get(sub(".db", "SYMBOL", annotation_package)), "all")
	}else if(missing(annot_env)){
		stop("Please give me a annotation!")
	}

	in.df <- read.delim(input, stringsAsFactors=FALSE)

	if(missing(probe_col)){
		message("Using probes in first column.")
		probe_col <- colnames(in.df)[1]
	}

	ann.df <- annotate(in.df, annot_env=annot_env, probe_col=probe_col, multiple=multiple)

	write.table(ann.df, output, sep="\t", row.names=FALSE, quote=FALSE)
}

if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	parameters <- arg

	if("annotation_package" %in% names(arg)){
		if(!suppressMessages(require(parameters$annotation_package, character.only=TRUE))) stop("Annotation package does not exists.")
	}else{
		if(!file.exists(arg$annotation_file)) stop("No such file: ", arg$annotation_file)
		annot_env <- get_annot_env(arg$annotation_file)
		parameters$annot_env <- annot_env
	}


	## RUN

	## Bibliotecas
	suppressMessages(library("GEOquery"))
	suppressMessages(library("jsonlite"))

	do.call(annotate_main, parameters)
}

