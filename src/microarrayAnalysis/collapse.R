#!/usr/bin/env Rscript

"Collapse

Usage: collapse.R INPUT OUTPUT --by-col=<value> [--method=<value>] [--annotation-cols=<value>...]

Input:
  INPUT                      the input file containing expression values with samples in columns and genes/probes in rows
  OUTPUT                     output file

Output:

Options:
  -h --help                  show this help message
  --version                  show program version
  --annotation-cols=<value>  annotation columns [default: ProbeName]
  --method=<value>           method (maxmean, minmean, median, colMeans)
  --by-col=<value>           column to group by

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
	return(in.df)
}

maxmean <- function(x){
	x[order(rowMeans(x), decreasing=TRUE)[1], ]
}

minmean <- function(x){
	x[order(rowMeans(x), decreasing=FALSE)[1], ]
}

colMedian <- function(x){
	apply(x, 2, median)
}

collapse <- function(in.df, method, by_col, drop=T){
	if(drop){
		tmp.df <- drop.cols(in.df, by_col)
	} else {
		tmp.df <- in.df
	}
	if(length(by_col) > 1){ 
		by_vals <- apply(in.df[, by_col], 1, paste, collapse="_")
	} else {
		 by_vals <- in.df[, by_col]
	}
	res <- tapply(1:nrow(in.df), by_vals, function(x){
		sub.df <- tmp.df[x, ]
		method(sub.df)
	})	
	return(do.call(rbind, res))
}

main_collapse <- function(input, output, method=c("maxmean", "minmean", "colMeans", "colMedian"), annotation_cols=c("ProbeName"), by_col){
	in.df <- read.delim(input)
	maintain.col <- which(annotation_cols %in% by_col)
	if(length(maintain.col) > 0){
		annotation_cols <- annotation_cols[-maintain.col]
	}
	exp.df <- collapse(drop.cols(in.df, annotation_cols), method=get(method), by_col=by_col)
	exp.df <- cbind(rownames(exp.df), exp.df)
	colnames(exp.df)[1] <- by_col
	
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

	do.call(main_collapse, parameters)
}

