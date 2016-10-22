#!/usr/bin/env Rscript

"Calculate correlation between lncRNAs and conding genes

Usage: cor_lnc.R --ovlp=<file> --exp=<file> --output=<file> --method=<value> --lnc-col=<value> --coding-col=<value> --symbol-col=<value> --probe-col=<value> --correction=<method> [--annotation-cols=<value>...]

Options:
  -h --help                   show this help message
  --version                   show program version
  --ovlp=<file>               file with overlaps between lncRNAs and coding
  --exp=<file>                expression file name
  --method=<value>            Correlation method (pearson or spearman) [default: spearman]
  --lnc-col=<value>           column on ovlp containing the lncRNA name [default: lnc_id]
  --coding-col=<value>        column on ovlp containing the coding [default: mrn_id]
  --symbol-col=<value>        names should match with lnc_col and coding_col [default: Symbol]
  --probe-col=<value>         column with probe names [default: ProbeName]
  --correction=<method>       method of correction p-value (holm, hochberg, hommel, bonferroni, BH, BY, fdr, none) [default: fdr]
  --output=<file>             output file name
  --annotation-cols=<value>   annotation columns

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc





## http://stats.stackexchange.com/questions/120199/calculate-p-value-for-the-correlation-coefficient
calc.p <- function(r, n){ 
    r <- abs(r) 
    r.1 <- matrix(mapply(function(x, y) {isTRUE(all.equal(x, y))}, r, 1),  
           nrow=nrow(r), ncol=ncol(r)) 
    r[which(r.1, arr.ind=T)] <- 1 
    t.val <- ( r * sqrt(n-2) ) / sqrt(1 - r**2) 
    p <- pt(t.val, df=n-2, lower.tail=FALSE)*2 
    return(p) 
}

p2r <- function(p, n){
    dfreedom <- n-2
    t.val <- qt(p/2, df=dfreedom, lower.tail=F)
    return(t.val/sqrt(t.val**2+dfreedom))
}



calc_correlation <- function(ovlp, exp_df, lnc_col, coding_col, symbol_col, probe_col, correction, annotation_cols){
	if(missing(annotation_cols)){
		annotation_cols <- c()
	}
	annotation_cols <- c(symbol_col, probe_col, annotation_cols)

	# indexes of coding and lnc in exp.df
	# (WARNING: row's order in exp_df shouldn't be modified)
	coding_rows <- which(exp_df[[symbol_col]] %in% ovlp[[coding_col]])
	lnc_rows <- which(exp_df[[symbol_col]] %in% ovlp[[lnc_col]])

	# convert exp_df to matrix
	exp_mat <- as.matrix(exp_df[, setdiff(colnames(exp_df), annotation_cols)])
	rownames(exp_mat) <- exp_df[[probe_col]]

	exp_trans <- t(exp_mat)

	message("Calculating correlation ...")
	cors <- cor(exp_trans[, coding_rows], exp_trans[, lnc_rows], method=method)

	# add probes related to lnc
	message("Looking for probes related to lncRNAs ...")
	ovlp <- merge(ovlp, exp_df[, annotation_cols], by.x=lnc_col, by.y=symbol_col, all.x=T)
	probe_col_lnc <- paste0(probe_col, "_lnc")
	colnames(ovlp)[grep(probe_col, colnames(ovlp))] <- probe_col_lnc

	# add probes related to coding genes
	message("Looking for probes related to coding genes ...")
	ovlp <- merge(ovlp, exp_df[, annotation_cols], by.x=coding_col, by.y=symbol_col, all.x=T)
	probe_col_coding <- paste0(probe_col, "_mrn")
	colnames(ovlp)[which(probe_col == colnames(ovlp))] <- probe_col_coding

	# remove rows without probes
	ovlp <- ovlp[-which(is.na(ovlp[, probe_col_coding]) | is.na(ovlp[, probe_col_lnc])), ]

	get_cor <- function(x, probe_coding=probe_col_coding, probe_lnc=probe_col_lnc, cor_mat=cors){
		coding <- x[probe_coding]
		lnc <- x[probe_lnc]
		return(cor_mat[coding, lnc])
	}

	message("Adding correlation to the data.frame ...")
	ovlp$correlation <- apply(ovlp, 1, get_cor)
	message("Calculating p-values ...")
	ovlp$p.value <- calc.p(matrix(ovlp$correlation), ncol(exp_mat))
	message("Adjusting p-values ...")
	ovlp$adj.p.value <- p.adjust(ovlp$p.value, method=correction)

	return(ovlp)
}


if (!interactive()) {
	# Get and check arguments.
	suppressMessages(library("docopt"))

	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	parameters <- arg

	## RUN
	suppressMessages(library("readr"))

	# parametros
	method <- arg[["method"]]  # Correlation method (pearson or spearman)
	ovlp_fname <- arg[["ovlp"]]  # file with overlaps between lncRNAs and coding
	exp_fname <- arg[["exp"]]  # expression file
	lnc_col <- arg[["lnc_col"]]  # column on ovlp containing the lncRNA name
	coding_col <- arg[["coding_col"]]  # column on ovlp containing the coding
	symbol_col <- arg[["symbol_col"]]  # names should match with lnc_col and coding_col
	probe_col <- arg[["probe_col"]]
	correction <- arg[["correction"]]

	if("annotation_cols" %in% names(arg)){
		annotation_cols <- arg[["annotation_cols"]]
	} else {
		annotation_cols <- c()
	}

	message("Reading overlap file ...")
	ovlp <- read_tsv(ovlp_fname)
	message("Reading expression file ...")
	exp_df <- read_tsv(exp_fname)

	corr <- calc_correlation(ovlp, exp_df, lnc_col, coding_col, symbol_col, probe_col, correction, annotation_cols)

	message("Writing output ...")
	write_tsv(corr, arg$output)
	message("DONE!")
}
