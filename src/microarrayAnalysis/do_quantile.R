#!/usr/bin/env Rscript

"Normalize a expression file

Usage: do_quantile.R INPUT OUTPUT  [--annotation-cols=<value>...] [--avoid-0] [--offset=<val>]

Input:
  INPUT                      the input file containing expression values with samples in columns and genes/probes in rows
  OUTPUT                     output file

Output:

Options:
  -h --help                  show this help message
  --version                  show program version
  --annotation-cols=<value>   annotation columns [default: ProbeName]
  --avoid-0                  if this flag is used and data is not log transformed I will sum abs(min) to all values 
  --offset=<val>             if avoid-0 is TRUE I will add a offset to all values

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

inner_quantile <- function(exp.df, method="quantile", annotation_cols=c("ProbeName"), avoid_0, offset_val){
	exp.mat <- normalizeBetweenArrays(as.matrix(drop.cols(exp.df, annotation_cols)), method=method)
	if(!is.log(exp.mat)){
		message("Data are not in log2 scale, I will apply log2 function.")
		min_val <- min(exp.mat)
		if(avoid_0 && min_val < 0){
			total_offset <- abs(min_val) + offset_val
			message("I will add ", round(total_offset, digits=2), " to all values before log transform.")
			exp.mat <- exp.mat + total_offset
		}
		exp.mat <- log2(exp.mat)
	}
	exp.df <- cbind(exp.df[, annotation_cols], as.data.frame(exp.mat))
	colnames(exp.df)[1:length(annotation_cols)] <- annotation_cols
	return(exp.df)
}

do_quantile <- function(input, output, method="quantile", annotation_cols=c("ProbeName"), avoid_0, offset_val, ...){
	exp.df <- read.delim(input)
	exp.df <- inner_quantile(exp.df, method, annotation_cols, avoid_0, offset_val)
	write.table(exp.df, output, quote=FALSE, sep="\t", row.names=FALSE)
}


if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	if(! "offset" %in% names(arg)){
		message("Setting offset to 0.")
		arg[["offset"]] <- 0
	}
	arg[["offset"]] <- as.numeric(arg[["offset"]])

	parameters <- arg

	## RUN
	library("limma")

	do.call(do_quantile, parameters)
}

