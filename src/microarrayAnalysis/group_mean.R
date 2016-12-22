#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

"Group mean

Usage: group_mean.R --expr-file=<file> --group-file=<file> --output=<file> [--group-col=<value>] [--samplename-col=<value>] 

Options:
  -h --help                   show this help message
  --version                   show program version
  --output=<file>             file to write the results
  --expr-file=<dir>           file to get expression values
  --group-file=<file>         file to get groups
  --samplename-col=<value>    column in group-file specifying the sample names [default: SampleName]
  --group-col=<value>         columns in group-file specifying the group [default: Group]

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
		return(in.df[, -rm.cols, drop=FALSE])
	}
	return(in.df)
}

if (!interactive()) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	## parameters
	out.file <- arg$output
	
	# read expression file
	exprs <- read.delim(file.path(arg$expr_file), stringsAsFactors=FALSE)
	# read the file with group information
	group.file <- read.delim(arg$group_file, stringsAsFactors=FALSE)
	# create a list of samples by groups
	groups <- split(group.file[, arg$samplename_col], group.file[, arg$group_col])
	# calculate the mean expression of each group and row
	groupmeans <- lapply(groups, function(g) apply(exprs[, g, drop=FALSE], 1, mean))
	# create a data.frame
	out <- as.data.frame(do.call(cbind, groupmeans))
	# add annotation columns in the output data.frame
	out <- cbind(drop.cols(exprs, group.file[, arg$samplename_col]), out)
	# write the table
	write.table(out, file.path(out.file), quote=FALSE, sep="\t", row.names=FALSE)
}

