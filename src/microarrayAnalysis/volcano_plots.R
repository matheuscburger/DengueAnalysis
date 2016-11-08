#!/usr/bin/env Rscript

"Volcano plots

Usage: volcano_plots.R --input=<file> --output=<file> [--log2fc-col=<value> --lfc-cut=<value> --pv-col=<value> --pv-cut=<value> --qv-col=<value> --qv-cut=<value> --p-chosen=<value> --cut-chosen=<value> (--y-min=<value> --y-max=<value>) (--x-min=<value> --x-max=<value>)]


Options:
  -h --help              show this help message
  --version              show program version
  --input=<file>         statistics file name
  --output=<file>        output filename (PDF)
  --log2fc-col=<value>   column containing log2(FC) [default: Log2FC]
  --lfc-cut=<value>      log2(FC) cutoff [default: 0]
  --pv-col=<value>       column containing p-value [default: LIMMA.rawp]
  --pv-cut=<value>       p-value cutoff [default: 1]
  --qv-col=<value>       column containing adjusted p-value [default: LIMMA.adjp]
  --qv-cut=<value>       adjusted p-value cutoff [default: 0.05]
  --p-chosen=<value>     p-value column to show [default: qv_col]
  --cut-chosen=<value>   cutoff on y-axis [default: qv_cut]
  --y-min=<value>        minimum value on y-axis
  --y-max=<value>        maximum value on y-axis
  --x-min=<value>        minimum value on x-axis
  --x-max=<value>        maximum value on x-axis


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

# parameters
volc_plot <- function(stats, log2fc_col, lfc_cut, pv_col, pv_cut, 
					  qv_col, qv_cut, out_fname, p_chosen, cut_chosen,
					  y_lim, x_lim){

	lfc.rows <- abs(stats[[log2fc_col]]) >= lfc_cut
	pv.rows <- stats[[pv_col]] <= pv_cut
	qv.rows <- stats[[qv_col]] <= qv_cut

	de.rows <- lfc.rows & pv.rows & qv.rows

	up.rows <- de.rows & stats[[log2fc_col]] > 0 
	down.rows <- de.rows & stats[[log2fc_col]] < 0 

	stats[, "Direction"] <- "Not Sig."
	stats[which(up.rows), "Direction"] <- "Up-regulated"
	stats[which(down.rows), "Direction"] <- "Down-regulated"

	volcano.colors <- c("Not Sig."="gray", 
						"Down-regulated"="#005417", 
						"Up-regulated"="#7D0002")

	pl <- stats %>% group_by(title) %>%
		do(
		   plots=ggplot(., aes_string(x=log2fc_col, y=paste0("-log10(",p_chosen,")"))) + 
			   geom_point(aes(color=Direction), alpha=0.5) +
				   geom_vline(xintercept=lfc_cut, color="#616161", linetype="dashed") +
				   geom_vline(xintercept=-lfc_cut, color="#616161", linetype="dashed") +
				   geom_hline(yintercept=-log10(cut_chosen), color="#616161", linetype="dashed") + theme_minimal() +
				   scale_colour_manual(values=volcano.colors) +
				   labs(x=expression("log2(Test/Control)"), y=expression("-log10(Significance)"), title=unique(.$title)) +
				   geom_text(data=. %>% count(Direction) %>%
							 filter(Direction != "Not Sig.") %>%
							 mutate(x=ifelse(Direction=="Down-regulated", -Inf, Inf)),
						 aes(label=n, x=x, color=Direction), y=Inf, vjust=2, hjust="inward", size=6, show.legend=F) +
		   			theme(legend.position="top")
				   )
	if(!missing(x_lim)){
		pl[["plots"]] <- lapply(pl[["plots"]], `+` , scale_x_continuous(limits=x_lim, oob=squish))
	}
	if(!missing(y_lim)){
		pl[["plots"]] <- lapply(pl[["plots"]], `+` , scale_y_continuous(limits=y_lim, oob=squish))
	}
	return(pl)
}


if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	par_names <- names(formals(volc_plot))
	parameters <- arg[which(names(arg) %in% par_names)]
	parameters[["lfc_cut"]] <- as.numeric(parameters[["lfc_cut"]])
	parameters[["pv_cut"]] <- as.numeric(parameters[["pv_cut"]])
	parameters[["qv_cut"]] <- as.numeric(parameters[["qv_cut"]])

	if(parameters[["p_chosen"]] == "qv_col"){
		parameters[["p_chosen"]] <- parameters[["qv_col"]]
	} else if(parameters[["p_chosen"]] == "pv_col"){
		parameters[["p_chosen"]] <- parameters[["pv_col"]]
	}

	if(parameters[["cut_chosen"]] == "qv_cut"){
		parameters[["cut_chosen"]] <- parameters[["qv_cut"]]
	} else if(parameters[["cut_chosen"]] == "pv_cut"){
		parameters[["cut_chosen"]] <- parameters[["pv_cut"]]
	} else {
		parameters[["cut_chosen"]] <- as.numeric(parameters[["cut_chosen"]])
	}

	if("y_min" %in% names(arg)){
		parameters[["y_lim"]] <- as.numeric(c(arg[["y_min"]], arg[["y_max"]]))
	}
	if("x_min" %in% names(arg)){
		parameters[["x_lim"]] <- as.numeric(c(arg[["x_min"]], arg[["x_max"]]))
	}

	## RUN
	parameters$stats <- read_tsv(arg$input)
	pl <- do.call(volc_plot, parameters)
	pdf(arg$output)
	lapply(pl[["plots"]], print)
	dev.off()
}
