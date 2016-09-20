#!/usr/bin/env Rscript

"Perform quality assessment

Usage: aqm.R (--input-exp=<value>|--input-dir=<value>) INPUT_ANNOT OUTPUT_DIR [--int-cols=<value>...] [--annotation-cols=<value>...]

Input:
  INPUT_ANNOT                the input file containing sample annotation
  OUTPUT                     output directory

Output:

Options:
  -h --help                  show this help message
  --version                  show program version
  --int-cols=<value>         Columns to show on aQM report
  --input-exp=<value>        the input file containing expression values with samples in columns and genes/probes in rows
  --input-dir=<value>        the input directory containing cel files
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
	return(in.df)
}

get_eset <- function(input_exp, input_annot, annotation_cols=c("ProbeName"), sample_name_col="Sample_geo_accession", int_cols=character(), ...){
	exp.df <- drop.cols(read.delim(input_exp), annotation_cols)
	sample.annot.tmp <- read.delim(input_annot)
	sel.cols <- c(sample_name_col, int_cols)
	sample.annot <- sample.annot.tmp[, sel.cols, drop=FALSE]
	na.cols <- which(apply(sample.annot, 2, function(x) sum(is.na(x))==length(x)))
	if(length(na.cols) > 0){
		sample.annot <- sample.annot[, -na.cols]
	}
	rownames(sample.annot) <- as.character(sample.annot[, sample_name_col])
	# remove samples from exp.df that are not present in sample annotation
	exp.df <- exp.df[, as.character(sample.annot[, sample_name_col])]
	sample.annot.obj <- AnnotatedDataFrame(sample.annot[colnames(exp.df), , drop=FALSE])
	message("Samples in expression file: ", paste(colnames(exp.df), sep=", ", collapse=", "))
	message("Samples in annotation file: ", paste(rownames(sample.annot.obj), sep=", ", collapse=", "))
	eset <- ExpressionSet(assayData=as.matrix(exp.df), phenoData=sample.annot.obj)
	return(eset)
}


get_affy <- function(input_dir, input_annot, output_dir, sample_name_col="Sample_geo_accession", cel_col="Sample_supplementary_file", int_cols=character(), ...){
	sample.annot.tmp <- read.delim(input_annot)
	files <- grep("\\.CEL", list.files(input_dir), ignore.case=T, value=T)
	sample.names <- sample.annot.tmp[, sample_name_col]

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

	sel.cols <- c(sample_name_col, int_cols)
	sample.annot <- sample.annot.tmp[, sel.cols, drop=FALSE]
	na.cols <- which(apply(sample.annot, 2, function(x) sum(is.na(x))==length(x)))
	if(length(na.cols) > 0){
		sample.annot <- sample.annot[, -na.cols]
	}
	rownames(sample.annot) <- as.character(sample.annot[, sample_name_col])
	sample.annot.obj <- AnnotatedDataFrame(sample.annot)
	message("Cel files: ", paste(cel.files, collapse=", "))
	message("Sample Names:  ", rownames(sample.annot.obj))
	affy.obj <- ReadAffy(celfile.path=input_dir, phenoData=sample.annot.obj, filenames=cel.files)
	return(affy.obj)
}

rm.tags <- function(x){
	return(gsub("<.*?>", "", x))
}

aqm <- function(input_exp, input_dir, input_annot, output_dir, annotation_cols=c("ProbeName"), sample_name_col="Sample_geo_accession", int_cols=character(0), ...){
	arguments <- as.list(match.call())[-1] # get args
	if(!missing(input_exp)){
		message("Get eset object...")
		eset <- do.call(get_eset, arguments)
		exp.obj <- eset
		int_cols <- int_cols[int_cols %in% varLabels(exp.obj)]
	}else if(!missing(input_dir)){
		message("Get AffyBatch object ...")
		affy.obj <- do.call(get_affy, arguments)
		exp.obj <- affy.obj
		int_cols <- int_cols[int_cols %in% varLabels(exp.obj)]
		prot <- pData(protocolData(exp.obj))
		if("ScanDate" %in% colnames(prot)){
			dates <- tryCatch(as.Date(prot[, "ScanDate"]),
							  error=function(x){as.Date(prot[, "ScanDate"], format="%m/%d/%Y")})
			if(length(unique(dates)) > 8){
				GroupByDate <- paste("day", cutree(hclust(dist(dates)), k=8), sep="_")
			}else{
				GroupByDate <- as.character(dates)
			}
			pData(exp.obj)[, "GroupByDate"] <- GroupByDate
			int_cols <- c("GroupByDate", int_cols)
		}
	}else{
		stop("Please specify input_exp or input_dir")
	}

	aqm.res <- arrayQualityMetrics(exp.obj, outdir=output_dir, intgroup=int_cols, force=T)

	aqm.tbl <- aqm.res$arrayTable
	colnames(aqm.tbl) <- rm.tags(colnames(aqm.tbl))
	write.table(aqm.tbl, file.path(output_dir, "outtable.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

	counts <- apply(aqm.tbl[, grep("\\*", colnames(aqm.tbl), value=T)], 1, function(x) sum(x=="x"))
	names(counts) <- aqm.tbl[, "sampleNames"]
	counts <- counts[counts > 0]
	writeLines(toJSON(list(total=length(grep("\\*", colnames(aqm.tbl))), counts=as.list(counts)), flatten=T, auto_unbox=T, pretty=T), file.path(output_dir, "outliers.json") )
}



if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library("docopt"))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	parameters <- arg

	## RUN
	suppressMessages(library("arrayQualityMetrics"))
	suppressMessages(library("Biobase"))
	suppressMessages(library("affy"))
	suppressMessages(library("jsonlite"))

	#save(file="aqm.RData", list=ls())
	do.call(aqm, parameters)
}

