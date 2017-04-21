#!/usr/bin/env Rscript

"Corrplot of fgsea output

Usage: corrplot_fgsea.R INPUT OUTPUT [--nes-suffix=<val> --pval-suffix=<val> --pathway-col=<val> --pv-cut=<val> --lim=<val>]

Input:
  INPUT                      the input file coming from fgsea
  OUTPUT                     output file (PDF)

Output:

Options:
  -h --help                  show this help message
  --version                  show program version
  --nes-suffix=<val>         suffix present on the columns containing NES [default: _NES]
  --pval-suffix=<val>        suffix present on the columns containing p-values [default: _padj]
  --pathway-col=<val>        column containing pathway names [default: pathway]
  --pv-cut=<val>             p-value cutoff [default: 0.1]
  --lim=<val>                limits for color legend [default: 2]

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

suppressMessages(library("corrplot"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))


get_mat <- function(tib, suffix, pathway_col){
	minus_pathway_col <- paste0("-", pathway_col) 
    res <- tib %>% 
	    dplyr::select_(minus_pathway_col) %>%
		dplyr::select(ends_with(suffix)) %>%
		as.matrix
	rownames(res) <- tib[[pathway_col]]
	colnames(res) <- sub(suffix, "", colnames(res))
	return(res)
}

print_corrplot <- function(gsea_results, corr_colors, nes_suffix, pval_suffix, pathway_col, pv_cut, out_fname){

	nes <- get_mat(gsea_results, nes_suffix, pathway_col)
	pvs <- get_mat(gsea_results, pval_suffix, pathway_col)

	# change row and col order
	pvs <- pvs[rownames(nes), colnames(nes)]

	# remove not significant rows
	rows_to_remove <- which(rowSums(pvs < pv_cut) == 0)

	pvs <- pvs[-rows_to_remove, ]
	nes <- nes[-rows_to_remove, ]

	# set to 0 the nes if p is not significant
	nes[which(pvs > pv_cut, arr.ind = T)] <- 0

	# change row order
	nes <- nes[order(rowSums(nes), decreasing=T), ]

	# define height
	# 52 modules for a height=7 is a good ratio
	height <- max(7, nrow(nes)*7/52) + 1

	# define cl.ratio
	clr <- (-0.07*height + 4.89)/44


	pdf(out_fname, height=height)
	corrplot(nes, 
			 col=corr_colors, is.corr=FALSE, addgrid.col="white", insig="blank",
			 pch.cex=0.5, pch.col="black", tl.col="black", tl.cex=0.5, cl.cex=0.4, cl.ratio=clr,
			 cl.pos="b", cl.align.text="l", mar=c(0,0,0,0), cl.lim=c(-4,4),
			 )
	dev.off()
}

if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))


	print(arg)
	## RUN

	#in_fname <- "results/fgsea/Log2FC/ReactomePathways_fixed/GSE13052.tsv"
	#out_fname <- "figures/CEMiTool_joined_corrplot/GSE13052.pdf" 
	#nes_suffix <- "_NES"
	#pval_suffix <- "_padj"
	#pathway_col <- "pathway"
	#pv_cut <- 0.1

	gsea_results <- read_tsv(arg$input)
	corr_colors <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
					"#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

	print_corrplot(gsea_results, corr_colors, arg$nes_suffix, arg$pval_suffix, arg$pathway_col, as.numeric(arg$pv_cut), arg$output)
}
