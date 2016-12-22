#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

"Z score

Usage: scale.R --expr-file=<file> --output=<file> [--annotation-cols=<value>...] 

Options:
  -h --help                   show this help message
  --version                   show program version
  --output=<file>             file to write the results
  --expr-file=<dir>           file to get expression values
  --annotation-cols=<value>   annotation columns [default: ProbeName]

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc



if (!interactive()) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	## parameters
	out.file <- arg$output
	ann.cols <- arg$annotation_cols

	exprs <- read.delim(file.path(arg$expr_file), stringsAsFactors=FALSE)

	columns <- colnames(exprs)
	value.columns <- columns[!columns %in% ann.cols]

	exprs[, value.columns] <- t(apply(exprs[, value.columns], 1, scale))

	write.table(exprs, file.path(out.file), quote=FALSE, sep="\t", row.names=FALSE)
}

