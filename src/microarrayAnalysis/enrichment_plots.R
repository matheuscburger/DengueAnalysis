#!/usr/bin/env Rscript

"Enrichment plots

Usage: enrichment_plots.R --input=<file> --output=<file> [--term-col=<value> --pv-col=<value> --rgb-color=<value> --pv-cut=<value> --max-length=<value> --max-show=<value> (--database-col=<value>|--title=<value>)]

Options:
  -h --help               show this help message
  --version               show program version
  --input=<file>          statistics file name
  --output=<file>         output filename (PDF)
  --term-col=<value>      column containing the term name [default: Term]
  --pv-col=<value>        columns containing the p-value [defaul: 'Adjusted P-value']
  --rgb-color=<value>     RGB color used in bars [default: #001773]
  --pv-cut=<value>        p-value threshold [default: 0.01]
  --max-length=<value>    maximum length of a term name [default: 60]
  --max-show=<value>      maximum number of rows to show [default: 10]
  --database-col=<value>  column containing the name of database [default: GeneSet]
  --title=<value>         a title

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

suppressMessages(library("dplyr"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggthemes"))
suppressMessages(library("scales"))
suppressMessages(library("stringr"))


my.squish <- function(...){
	return(squish(..., only.finite=FALSE))
}

plot_ora <- function(es, name, order_by, plot_title, graphcolor, pv_cut, maxlength, max_show, wrap_width=30) {

	if(max_show > nrow(es)){
		max_show <- nrow(es)
	}
	
	es <- head(es[order(es[[order_by]]), ], n=max_show)

	# limits name length
	rows_idx <- which(nchar(es[[name]]) > maxlength) # rows exceding maxlength
	es[[name]][rows_idx] <- paste0(strtrim(gsub("_", " ", es[[name]][rows_idx]), maxlength), "...")
	es[[name]] <- str_wrap(es[[name]], width = wrap_width)

	# order bars
	lvls <- es[[name]][order(-log10(es[[order_by]]), decreasing=F)]
	es[[name]] <- factor(es[[name]], levels=lvls)

	es[["alpha"]] <- 1
	es[["alpha"]][es[, order_by] > pv_cut] <- 0

	# Avoid 0's
	es[[order_by]][es[[order_by]] > 0.8] <- 0.8

	order_by.name <- paste0("-log10(", order_by, ")")

	# plot
	pl <- ggplot(es, aes_string(x=name, y=paste0('-log10(`', order_by, '`)'), alpha="alpha", fill=paste0('-log10(`', order_by, '`)'))) + 
		geom_bar(stat="identity") + 
		theme(axis.text=element_text(size=8), legend.title=element_blank()) +
		coord_flip() + scale_alpha(range=c(0.4, 1), guide="none") +
		labs(y=order_by.name, fill=order_by.name, title=plot_title, x="") +
		geom_hline(yintercept=-log10(pv_cut), colour="grey", linetype="longdash") + 
		scale_fill_gradient(low="gray", high=graphcolor, limits=c(2, 5), oob=my.squish) +
		theme_minimal()
	return(pl)
}



if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library("docopt"))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))
	if(is.null(arg$pv_col)){
		arg$pv_col <- "Adjusted P-value"
	}

	## RUN
	in_df <- read_tsv(arg$input)
	if(!is.null(arg$title)){
		in_df[, "GeneSet"] <- arg$title
		arg$database_col <- "GeneSet"
	} else {
		in_df[, "GeneSet"] <- ""
		arg$database_col <- "GeneSet"
    } 

	pl <- in_df %>% group_by_(arg$database_col) %>% 
		do(
			plots=plot_ora(., name=arg$term_col, order_by=arg$pv_col, 
							  plot_title=unique(.[[arg$database_col]]),
							  graphcolor=arg$rgb_color, pv_cut=as.numeric(arg$pv_cut), 
							  maxlength=as.numeric(arg$max_length), max_show=as.numeric(arg$max_show),
							  wrap_width=30)
		)

	pdf(arg$output)
	lapply(pl[["plots"]], print)
	dev.off()
}
