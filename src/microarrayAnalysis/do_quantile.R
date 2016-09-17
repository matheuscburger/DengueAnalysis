#!/usr/bin/env Rscript

"Normalize a expression file

Usage: do_quantile.R INPUT OUTPUT  [--annotation-cols=<value>...]

Input:
  INPUT                      the input file containing expression values with samples in columns and genes/probes in rows
  OUTPUT                     output file

Output:

Options:
  -h --help                  show this help message
  --version                  show program version
  --annotation-cols=<value>   annotation columns [default: ProbeName]

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

drop.cols <- function(in.df, cols){
	rm.cols <- which(colnames(in.df) %in% cols)
	if(length(rm.cols) > 0){
		return(in.df[, -rm.cols])
	}
	return(in.df[, -rm.cols])
}

is.log <- function(mat, limit=1000){
	return(!max(mat, na.rm=TRUE) > 1000)
}

inner_quantile <- function(exp.df, method="quantile", annotation_cols=c("ProbeName")){
	exp.mat <- normalizeBetweenArrays(as.matrix(drop.cols(exp.df, annotation_cols)), method=method)
	if(!is.log(exp.mat)){
		exp.mat <- log2(exp.mat)
	}
	exp.df <- cbind(exp.df[, annotation_cols], as.data.frame(exp.mat))
	colnames(exp.df)[1:length(annotation_cols)] <- annotation_cols
	return(exp.df)
}

do_quantile <- function(input, output, method="quantile", annotation_cols=c("ProbeName"), ...){
	exp.df <- read.delim(input)
	exp.df <- inner_quantile(exp.df, method, annotation_cols)
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
	library("limma")

	do.call(do_quantile, parameters)
}

