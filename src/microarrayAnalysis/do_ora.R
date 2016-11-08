#!/usr/bin/env Rscript

options(stringsAsCharacter=FALSE)

"Do ORA

Usage: do_ora.R OUTPUT (--stats=<file>|--genes=<file>) --gmt=<file> [--pv-col=<value> --padj-col=<value> --lfc-col=<value> --comp-col=<value> --gene-col=<value> --pv-cut=<value> --padj-cut=<value> --lfc-cut=<value>]


Output:
  OUTPUT                     output file

Options:
  -h --help                  show this help message
  --version                  show program version
  --stats=<file>             the input file containing the statistics (format: TSV)
  --genes=<file>             the input file containing the genes (format: a list )
  --gmt=<file>               GMT (gene matrix transposed) file
  --comp-col=<value>         column in INPUT containing the name of comparison [default: title] 
  --gene-col=<value>         column in INPUT containing the gene [default: Symbol]
  --pv-col=<value>           column in INPUT containing the raw p-value [default: LIMMA.rawp]
  --padj-col=<value>         column in INPUT containing the adjusted p-value [default: LIMMA.adjp]
  --lfc-col=<value>          column in INPUT containing the log2 of fold-change [default: Log2FC]
  --pv-cut=<value>           p-value cutoff [default: 0.05]
  --padj-cut=<value>         adjusted p-value cutoff [default: 0.05]
  --lfc-cut=<value>          log2(FC) cutoff [default: 1]

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc



# arquivo GMT
read.gmt <- function(fname){
    res <- list(genes=list(), desc=list())
    gmt <- file(fname)
    gmt.lines <- readLines(gmt)
    close(gmt)
    gmt.list <- lapply(gmt.lines, function(x) unlist(strsplit(x, split="\t")))
    gmt.names <- sapply(gmt.list, '[', 1)
    gmt.desc <- lapply(gmt.list, '[', 2)
    gmt.genes <- lapply(gmt.list, function(x){x[3:length(x)]})
    names(gmt.desc) <- names(gmt.genes) <- gmt.names
    return(list(genes=gmt.genes, desc=gmt.desc))
}


prepare.gmt <- function(gmt){
	res <- list()
	res[["term2gene"]] <- do.call(rbind, lapply(names(gmt[["genes"]]), 
												function(n) cbind(Term=n, Gene=gmt[["genes"]][[n]])))
	res[["term2name"]] <- do.call(rbind, lapply(names(gmt[["desc"]]), 
												function(n) cbind(Term=n, Name=gmt[["desc"]][[n]])))
	return(res)
}

do.ora.list <- function(topgenes, gmt.list, allgenes, ...){
	if(missing(allgenes)){
		message("Using all genes in GMT file as universe.")
		allgenes <- unique(gmt.list[["term2gene"]][, "Gene"])
	}
	enriched <- enricher(gene = topgenes,
						 pvalueCutoff = 1,
						 qvalueCutoff = 1,
						 universe = allgenes,
						 TERM2GENE = gmt.list[['term2gene']],
						 TERM2NAME = gmt.list[['term2name']])
	result <- cbind(title=NA, direction=NA, enriched@result)
	return(result)
}

do.ora.stats <- function(stat_df, gmt.list, pv_cut, padj_cut, lfc_cut, gene_col, comp_col, lfc_col, pv_col, padj_col, ...){


	comparisons <- stat_df %>% select(get(comp_col)) %>% unique %>% unlist
	all.comp <- data.frame()

	for(cmp in comparisons){
		message("Comparison: ", cmp)

        # filter for rows
		comp.rows <- stat_df[, comp_col] == cmp
		allgenes <- stat_df[comp.rows, gene_col] %>% unique %>% collect %>% .[[gene_col]]
		pv.rows <- stat_df[, pv_col] < pv_cut
		padj.rows <- stat_df[, padj_col] < padj_cut
		up.rows <- stat_df[, lfc_col] >= lfc_cut
		down.rows <- stat_df[, lfc_col] <= -(lfc_cut)
        # messages about the filtering
		message("Table with statistics: ")
		print(head(stat_df[comp.rows, ]))
		message("Number of rows after p-value filter: ", sum(pv.rows))
		message("Number of rows after Adjusted P-value filter: ", sum(padj.rows))

		for(direction in c("up", "down")){
            # Get topgenes
			message("Direction: ", direction)
			if(direction == "up"){
				topgenes <- stat_df[which(comp.rows & pv.rows & padj.rows & up.rows), gene_col] %>% 
				collect %>% .[[gene_col]]
			} else {
				topgenes <- stat_df[which(comp.rows & pv.rows & padj.rows & down.rows), gene_col] %>% 
				collect %>% .[[gene_col]]
			}
			topgenes <- unique(topgenes)
			topgenes <- topgenes[!is.na(topgenes)]

			if(length(topgenes) == 0){
				message("No significant gene.")
				next
			} else {
				message("Genes used (", length(topgenes),"): ", paste(topgenes, collapse=", "))
			}

			# ORA
			enriched <- enricher(gene = topgenes,
								 pvalueCutoff = 1,
								 qvalueCutoff = 1,
								 universe = allgenes,
								 TERM2GENE = gmt.list[['term2gene']],
								 TERM2NAME = gmt.list[['term2name']])
			result <- cbind(cmp, direction=direction, enriched@result)
			colnames(result)[1] <- comp_col
			all.comp <- rbind(all.comp, result)
		}
	}
	return(all.comp)
}


if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library("docopt"))
	suppressMessages(library("clusterProfiler"))
	suppressMessages(library("readr"))
	suppressMessages(library("tidyr"))
	suppressMessages(library("dplyr"))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

    print(arg)
	# convert to numeric
	for(n in c("pv_cut", "padj_cut", "lfc_cut")){
		arg[[n]] <- as.numeric(arg[[n]])
	}

	# GMT
	gmt <- read.gmt(arg$gmt)
	gmt.list <- prepare.gmt(gmt)

	# arguments
	parameters <- arg
	parameters[["gmt.list"]] <- gmt.list


	if ("stats" %in% names(arg)){
		do.ora <- do.ora.stats
		suppressMessages(stat_df <- read_tsv(arg$stats))
		parameters[["stat_df"]] <- stat_df
	} else {
		do.ora <- do.ora.list
		parameters[["topgenes"]] <- readLines(arg[["genes"]])
	}

	## RUN
	all.comp <- do.call(do.ora, parameters)
	write.table(all.comp, arg$output, row.names=FALSE, sep="\t", quote=FALSE)
}

