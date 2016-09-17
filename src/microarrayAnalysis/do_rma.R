#!/usr/bin/env Rscript

"Do RMA

Usage: do_rma.R --sup-dir=<dir> --rma-file=<dir> --sample-annot=<file> [--sample-name-col=<value>] [--cdf-package=<pack>]

Options:
  -h --help                    show this help message
  --version                    show program version
  --cdf-package=<pack>         cdf package
  --sup-dir=<dir>              directory to get CEL files
  --rma-file=<file>            file to write normalized RMA 
  --sample-annot=<file>        directory to write sample annotation file
  --sample-name-col=<value>    name of the column containing the sample names [default: Sample_geo_accession]

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

#' Do log-scale Robust Multi-Array Analysis
#' RMA uses a linear model of the form: T(PM_ij) = e_i + a_j + epsilon_ij
#' T is the transformation that background corrects, normalizes and log PM intensities
#' e_i is the log2 scale expression value found on array i
#' a_j is the log scale affinity effects for probe j 
#' epsilon_ij is the error for probe j in array i
#' A robust linear-fitting procedure is used to estimate log-scale expression values e_i
#' (Statistics and Data Analysis for Microarrays Using R and Bioconductor - Sorin Draghici)
#' 
do_rma <- function(cdf.pack, sup.dir, sample.annot, rma.file, sample.name.col, ...){

	# get all samples
	sample.annotation <- read.delim(sample.annot, stringsAsFactors=FALSE)	

	sample.names <- sample.annotation[, sample.name.col]

	cel.dir <- file.path(sup.dir)
	files <- grep("cel", list.files(cel.dir, full.names=T), value=T, ignore.case=TRUE)
	cel.files <- sapply(sample.names, function(x){
		res <- grep(x, files, value=T)
		if(length(res) != 1){
			return(NA)
		}else{
			return(res)
		}
	})
	if(any(is.na(cel.files))){
		notfound <- paste(sample.names[which(is.na(cel.files))], collapse=", ")
		stop("Some cel files not found! Samples:", notfound)
	}
	message("Cel files: ", paste(cel.files, collapse=", "))
	message("Sample names: ", paste(sample.names, collapse=", "))

	message("Doing RMA ...")
	if(missing(cdf.pack)){
		rma.eset <- justRMA(filenames=basename(cel.files), celfile.path=cel.dir, sampleNames=sample.names)
	}else{
		rma.eset <- justRMA(filenames=basename(cel.files), celfile.path=cel.dir, sampleNames=sample.names, cdfname=cdf.pack)
	}
	rma.df <- as.data.frame(exprs(rma.eset))
	rma.df <- cbind(ProbeName=rownames(rma.df), rma.df)
	message("Writing results ...")
	write.table(rma.df, file.path(rma.file), quote=FALSE, sep="\t", row.names=FALSE)
}


if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library("docopt"))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-|_', '.', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	parameters <- arg

	if("cdf.package" %in% names(arg)){
		cdf.pack <- arg$cdf.package
		parameters[['cdf.pack']] <- cdf.pack
	}


	## RUN

	## Bibliotecas
	suppressMessages(library("affy"))
	suppressMessages(library("jsonlite"))

	message("Calling do_rma")
	do.call(do_rma, parameters)
	message("do_rma finished")

}

