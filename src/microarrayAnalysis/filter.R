#!/usr/bin/env Rscript

"Filter

Usage: filter.R INPUT OUTPUT --method=<value> (--prop=<value>|--rows=<value>) [--annotation-cols=<value>...]

Input:
  INPUT                      the input file containing expression values with samples in columns and genes/probes in rows
  OUTPUT                     output file

Output:

Options:
  -h --help                  show this help message
  --version                  show program version
  --annotation-cols=<value>  annotation columns [default: ProbeName]
  --prop=<value>             proportion (a value between 0 and 1)
  --rows=<value>             number of rows to be mantained
  --method=<value>           method to be used for order rows (sd for standard deviation or mean for average)

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

filter.prop <- function(x, prop, fun=mean, annotation_cols=c("ProbeName")){
    rows <- floor(prop * nrow(x))
    return(filter.rows(x, rows, fun, annotation_cols))
}

filter.rows <- function(x, rows, fun=sd, annotation_cols=c("ProbeName")){
    rows <- min(rows, nrow(x))
	dropped <- drop.cols(x, annotation_cols)
	if(!all(apply(dropped, 2, is.numeric))){
		stop("Values are not numeric!")
	}
    val <- apply(dropped, 1, fun)
    sel.rows <- order(val, decreasing=TRUE)[1:rows]
    return(x[sel.rows, ])
}


drop.cols <- function(in.df, cols){
	rm.cols <- which(colnames(in.df) %in% cols)
	if(length(rm.cols) > 0){
		return(in.df[, -rm.cols])
	}
	return(in.df[, -rm.cols])
}


main_filter <- function(input, output, prop, rows, method, annotation_cols=c("ProbeName"), ...){
	if(!missing(prop) && !missing(rows)){
		stop("Please specify prop or rows, not both.")
	}

	if(!method %in% c("mean", "sd")){
		stop("Wrong method!")
	}
	in.df <- read.delim(input)

	if(!missing(prop)){
		exp.df <- filter.prop(in.df, as.numeric(prop), fun=get(method), annotation_cols)
	}else if(!missing(rows)){
		exp.df <- filter.rows(in.df, as.numeric(rows), fun=get(method), annotation_cols)
	}else{
		stop("Please, specify prop or rows.")
	}
	
	write.table(exp.df, output, quote=FALSE, sep="\t", row.names=FALSE)
}


if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	parameters <- arg

	## RUN

	do.call(main_filter, parameters)
}

