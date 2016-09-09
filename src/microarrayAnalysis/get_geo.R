#!/usr/bin/env Rscript

"Download files from geo

Usage: get_geo.R GEO_ID [PLAT_ID] [--author-dir=<dir> --sup-dir=<dir> --sample-annot-dir=<dir> --probe-annot-dir=<dir> --raw-dir=<dir>]

Input:
  GEO_ID                     a study ID from GEO database
  PLAT_ID                    a platform ID from GEO database

Output:

Options:
  -h --help                  show this help message
  --version                  show program version
  --annotation-file=file     annotation file containing the probe id in the first column e symbol gene in the second column
  --annotation-package=pack  annotation package (AnnotationDbi)
  --author-dir=dir           directory to write author normalized data
  --sup-dir=dir              directory to write supplemental files
  --raw-dir=dir              directory to write raw files
  --sample-annot-dir=dir     directory to write sample annotation file
  --probe-annot-dir=dir      directory to write probe annotation file


Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc


#' The command getGEO from GEOquery package try to divide some columns in phenoData. This function
#' will regroup this columns and add a new column with GEO id.
#'
#' @param pd	the data.frame returned by pData function (sample annotation).
#' @param geo_id	GEO id
#' @return a new data frame with some columns regrouped and GEO id
regroup_pd <- function(pd, geo_id){
	## Regroup columns (example: characteristics_ch1, characteristics_ch1.1, characteristics_ch1.2, ...)
	cols_to_cat <- unique(gsub("\\.\\d", "", grep("\\.\\d$", colnames(pd), value=T)))
	for(mycol in cols_to_cat){
		indices <- grep(mycol, colnames(pd))
		newcol <- apply(pd[, indices], 1, paste, collapse=";")
		pd <- pd[, -indices]
		pd[, mycol] <- newcol
		pd <- cbind(GSE=geo_id, pd)
	}
	return(pd)
}

#' Get data from GEO. 
#'
#' @param geo_id	GEO study d
#' @param plat_id	GEO platform id
#' @param author_dir	directory to write author normalized data
#' @param sup_dir	directory to write supplemental files
#' @param sample_annot_geo_dir	directory to write sample annotation files
#' @param probe_annot_dir	directory to write probe annotation files (from GEO)
get_geo <- function(geo_id, plat_id, author_dir, sup_dir, sample_annot_dir, probe_annot_dir, raw_dir, ...){

	message("Getting supplemental files ...")
	if(!missing(sup_dir)){
		supp_log <- getGEOSuppFiles(geo_id, baseDir=file.path(sup_dir))

		message("Uncompressing tar files ...")
		tarfiles <- grep("\\.tar$", rownames(supp_log), value=T, ignore.case=TRUE)
		for(tfile in tarfiles){
			res <- untar(tfile, exdir=file.path(sup_dir, geo_id))
		}
	}

	message("Getting ExpressionSet from GEO ...")
	geo <- getGEO(geo_id)
	for(eset in geo){
		if(!missing(plat_id)){
			if(plat_id != annotation(eset)){
				next
			}
		}else{
			plat_id <- annotation(eset)
		}

		name <- paste0(geo_id, "_", annotation(eset))
		fname <- paste0(name, ".tsv")

		# move raw files from suplementary files directory to raw files directory
		if(!missing(raw_dir)){
			dir.to.create <- file.path(raw_dir, paste(geo_id, plat_id, sep="_"))
			if(!dir.exists(dir.to.create)){
				flag <- dir.create(dir.to.create)
				if(!flag) stop("Error creating raw dir.")
			}
			supcols <- grep("supplementary_file", colnames(pData(eset)), value=TRUE) # columns with supplementar material
			supfiles <- basename(as.character(unlist(pData(eset)[, supcols])))

			allsup <- list.files(file.path(sup_dir, geo_id))

			missingfiles <- setdiff(supfiles, allsup)
			presentfiles <- intersect(supfiles, allsup)

			if(length(missingfiles) > 0){
				warning("Some supplementary files were not found! Files: ", paste(missingfiles, collapse=","))
			}
			flag <- file.rename(file.path(sup_dir, geo_id, presentfiles), file.path(dir.to.create, presentfiles))
			if(!all(flag)) stop("Error moving files.")
		}

		# author expression data
		if(!missing(author_dir)){
			exp_df <- as.data.frame(exprs(eset))
			exp_df <- cbind(ProbeName=rownames(exp_df), exp_df)
			write.table(exp_df, file.path(author_dir, fname), quote=FALSE, sep="\t", row.names=FALSE)
			message("Author files created.")
		}

		# sample annotation
		if(!missing(sample_annot_dir)){
			pd <- regroup_pd(pData(eset), geo_id)
			write.table(pd, file.path(sample_annot_dir, fname), quote=FALSE, sep="\t", row.names=FALSE)
			message("Sample annotation files created.")
		}

		# probe_annotation
		if(!missing(probe_annot_dir)){
			write.table(pData(featureData(eset)), file.path(probe_annot_dir, fname), quote=FALSE, sep="\t", row.names=FALSE)
			message("Probe annotation files created.")
		}

	}
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

	## Bibliotecas
	suppressMessages(library("GEOquery"))
	suppressMessages(library("jsonlite"))


	do.call(get_geo, parameters)
}

