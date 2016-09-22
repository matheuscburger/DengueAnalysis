#!/usr/bin/env Rscript

"Do comparisons

Usage: do_comparisons.R --comp-file=<file> --norm-file=<file> --sample-annot=<file> --stat=<file> [--sample-name-col=<value>] [--annotation-cols=<value>...]

Options:
  -h --help                   show this help message
  --version                   show program version
  --sample-annot=<dir>        file to get sample annotation file
  --stat=<dir>                file to write statistics
  --norm-file=<dir>           file to get normalized expression values
  --comp-file=<dir>           file describing all comparisons
  --annotation-cols=<value>   annotation columns [default: ProbeName]
  --sample-name-col=<value>   name of the column containing the sample names [default: Sample_geo_accession]

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

forbidden_chars <- function(x, accepted.chars=c('\\(', '\\)', '\\s+', '\'', '[A-Za-z0-9_.]', ',')){
	res <- x
	for(acc in accepted.chars){
		res <- str_replace_all(res, acc, "")
	} 
	return(nchar(res) > 0)
}

convert_search <- function(x){
	if(forbidden_chars(x)){
		stop("Forbidden character.")
	}
	res <- x
	res <- str_replace_all(res, "\\(([A-Za-z0-9_.]+),\\s+('?[A-Za-z0-9_.]+'?)\\)", "(\\1 == \\2)")
	res <- str_replace_all(res, "(\\s*OR\\s*)", " | ")
	res <- str_replace_all(res, "(\\s*AND\\s*)", " & ")
	return(res)
}

sel_lines <- function(in.df, search.exprs){
	search.str <- convert_search(search.exprs)
	return(which(eval(parse(text=search.str), envir=in.df)))
}

#' Calculate log2(Fold-change)
#' @param x a numeric vector with expression values
#' @param samp.group a list specifying test and control samples
#' 
#' @return log2(fold-change)
calc.lfc <- function(x, samp.group){
	lfc <- mean(as.numeric(x[samp.group[["test"]]]), na.rm=TRUE) - mean(as.numeric(x[samp.group[["control"]]]), na.rm=TRUE)
	return(lfc)
}

#' Get p-values
#' @param x a numeric vector with expression values
#' @param samp.group a list specifying test and control samples
#' @param paired is test paired ?
#' @param test function that calculate the test
#'
#' @return a p-value
calc.test <- function(x, samp.group, paired=FALSE, test=t.test){
	my.t <- try(test(x=as.numeric(x[samp.group[["control"]]]), y=as.numeric(x[samp.group[["test"]]]), paired=paired), silent=T)
	if (is(my.t, "try-error")) {
		res <- NA
	}else{
		res <- my.t$p.value
	}
	return(res)
}

#' Do comparisons
#' @param exprs a data.frame containing expression values (rows represents genes and columns represents samples)
#' @param sample.ann a data.frame containing information from each sample (rows)
#' @param control.name name of the group containing the control samples
#' @param test.name name of the group containing the test samples
#' @param title title of comparison
#' @param class.col class specifying groups
#'
#' @return a data.frame containing statistics
do.comp <- function(exprs, sample.ann, comp, annotation_cols=c("ProbeName", "Symbol"), sample_name_col){

	# nome da comparacao
	comp.name <- gsub(" ", ".", comp['title'])
	message("Comparison name: ", comp.name, "\n")

	# amostras de cada grupo
	control.samples <- sample.ann[sel_lines(sample.ann, comp['control']), sample_name_col]
	test.samples <- sample.ann[sel_lines(sample.ann, comp['test']), sample_name_col]
	samp.group <- list("test"=test.samples, "control"=control.samples)

	message("Control samples: ", paste(control.samples, collapse=", "), "\n")
	message("Test samples: ", paste(test.samples, collapse=", "), "\n")

	# subset exprs
	exprs <- exprs[, c(annotation_cols, control.samples, test.samples)]
	message("head(exprs) ...")
	sink(stderr())
	print(exprs[1:4, 1:4])
	sink()

	# data.frame com estatisticas
	stats <- exprs[, annotation_cols, drop=FALSE]

	# Calcula media
	message("Calculating average")
	stats[rownames(exprs), paste0("Average")] <- apply(exprs[, c(control.samples, test.samples)] , 1, mean, na.rm=TRUE)
	message("head(Average) ...")
	sink(stderr())
	print(head(stats[rownames(exprs), paste0("Average")]))
	sink()

	# Calcula fold-change
	message("Calculating fold-change")
	stats[rownames(exprs), paste0("Log2FC")] <- apply(exprs[, c(control.samples, test.samples)] , 1, calc.lfc, samp.group)
	message("head(Log2FC) ...")
	sink(stderr())
	print(head(stats[rownames(exprs), paste0("Log2FC")]))
	sink()

	# Calcula teste T
	message("Doing t-test")
	stats[rownames(exprs), paste0("Ttest.rawp")] <- apply(exprs[, c(control.samples, test.samples)], 1, calc.test, samp.group, test=t.test)
	stats[rownames(exprs), paste0("Ttest.adjp")] <- p.adjust(stats[rownames(exprs), paste0("Ttest.rawp")], method="fdr")
	message("head(Ttest.rawp) ...")
	sink(stderr())
	print(head(stats[rownames(exprs), paste0("Ttest.rawp")]))
	sink()

	# Calcula Wilcoxon
	message("Doing Wilcoxon ...")
	stats[rownames(exprs), paste0("Wilcoxon.rawp")] <- apply(exprs[, c(control.samples, test.samples)], 1, calc.test, samp.group, test=wilcox.test)
	stats[rownames(exprs), paste0("Wilcoxon.adjp")] <- p.adjust(stats[rownames(exprs), paste0("Wilcoxon.rawp")], method="fdr")
	message("head(Wilcoxon.rawp) ...")
	sink(stderr())
	print(head(stats[rownames(exprs), paste0("Wilcoxon.rawp")]))
	sink()

	## Calcula SAM
	#resp <- c(rep(1, length(control.samples)), rep(2, length(test.samples)))
	#samfit <- SAM(exprs[,c(control.samples, test.samples)], resp, resp.type="Two class unpaired", fdr.output=params$qv.cut, geneid=rownames(exprs), logged2=TRUE, nperm=1000)
	#stats[samfit[["siggenes.table"]][["genes.up"]][, "Gene Name"], paste0("SAM.fdr")] <- as.numeric(samfit[["siggenes.table"]][["genes.up"]][, "q-value(%)"])/100
	#stats[samfit[["siggenes.table"]][["genes.down"]][, "Gene Name"], paste0("SAM.fdr")] <- as.numeric(samfit[["siggenes.table"]][["genes.down"]][, "q-value(%)"])/100

	# Calcula LIMMA
	message("Doing LIMMA ... ")
	Class <- factor(c(rep("control", length(control.samples)), rep("test", length(test.samples))))
	design <- model.matrix(~0+Class)
	colnames(design) <- levels(Class)
	cm <- makeContrasts(TvsC = test-control, levels=design)
	colnames(cm) <- comp.name
	fit <- lmFit(exprs[,c(control.samples, test.samples)], design)
	fit2 <- contrasts.fit(fit, cm)
	fit2 <- eBayes(fit2)
	l.stats <- topTable(fit2, comp.name, number=Inf)	# limma stats
	stats[rownames(l.stats), paste0("LIMMA.rawp")] <- l.stats[, "P.Value"]
	stats[rownames(l.stats), paste0("LIMMA.adjp")] <- l.stats[, "adj.P.Val"]
	message("head(LIMMA.rawp) ...")
	sink(stderr())
	print(head(stats[rownames(exprs), paste0("LIMMA.rawp")]))
	sink()

	return(stats)
}

if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	# carrega pacotes
	suppressMessages(library("limma"))
	suppressMessages(library("jsonlite"))
	suppressMessages(library("samr"))
	suppressMessages(library("stringr"))

	## parameters
	comp.file <- arg$comp_file
	stat.file <- arg$stat
	annotation_cols <- arg$annotation_cols


	sample.ann.file <- arg$sample_annot
	bname <- sub(".tsv", "", basename(sample.ann.file))

	# arquivo de anotacao das amostas
	sample.ann <- read.delim(sample.ann.file, stringsAsFactor=FALSE)
	rownames(sample.ann) <- sample.ann[, arg$sample_name_col]

	exprs <- read.delim(file.path(arg$norm_file), stringsAsFactors=FALSE)

	# le arquivo com comparacoes
	comp <- read.delim(comp.file, stringsAsFactor=FALSE)
	sink(stderr())
	print(comp)
	sink()
	comp.res <- list()
	for( i in 1:nrow(comp) ){
		stat.df <- do.comp(exprs, sample.ann, comp[i,], annotation_cols, sample_name_col=arg$sample_name_col)
		stat.df <- cbind(title=comp[i, "title"], stat.df)
		comp.res[[comp[i, "title"]]] <- stat.df
	}
	all.stats <- do.call(rbind, comp.res)
	write.table(all.stats, file.path(stat.file), quote=FALSE, sep="\t", row.names=FALSE)
}

