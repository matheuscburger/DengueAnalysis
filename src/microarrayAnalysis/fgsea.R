#!/usr/bin/env Rscript

"FGSEA of mean Z-score

Usage: fgsea.R --expression=<file> --es=<file> --nes=<file> --pval=<file> --gmt=<file> --sample-annotation=<file> --annotation-cols=<value> --symbols=<value> --sample-name-col=<value> --class-col=<value> [--annotation-cols=<value>...]

Options:
  -h --help                  show this help message
  --version                  show program version
  --expression=<file>        expression file
  --es=<file>                output file containing enrichment scores
  --nes=<file>               output file containing normalized enrichment scores
  --pval=<file>              output file containing p-values
  --gmt=<file>               gmt file
  --sample-annotation=<file> sample annotation file name
  --annotation-cols=<value>  annotation columns
  --symbols=<value>          column containing symbols
  --sample-name-col=<value>  columns containing sample name in sample annotation file
  --class-col=<value>        column containing class in sample annotation file

Authors:
  Matheus Carvalho Burger - burger at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

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
    res[['genes']] <- gmt.genes
    res[['desc']] <- gmt.desc
    return(res)
}


doFastGSEA <- function(exp.gsea, template.df, GS, ranks = F, class_col){ # from CEMiTool with modifications
    
	message("Running function doFastGSEA ...")
    Temp <- template.df[which(!is.na(template.df[[class_col]])), ] # removing rows with class == NA
    Temp[[class_col]] <- as.character(Temp[[class_col]])
    classes <- as.vector(unique(Temp[, class_col]))

	exp.gsea <- exp.gsea[, rownames(Temp)]
    
    Zexp.gsea <- data.frame(t(scale(t(exp.gsea), center=TRUE, scale=TRUE))) # transforma em tabela d Z-scores
    
    gseaList <- list()
	
    
    for(j in 1:length(classes)){
        curr_class <- classes[j]
        
        class_samples <- rownames(subset(Temp, subset=Temp[[class_col]]==curr_class))
        
		message("Class ", j, " ", curr_class, " ...")
        message("Samples: ", paste(class_samples, collapse=", "))
        
        if(ranks){
            geneList <- rank(apply(exp.gsea[, class_samples], 1, mean))
            geneList <- sort(geneList, decreasing = T)
        }else{
            geneList <- apply(Zexp.gsea[, class_samples], 1, mean)
            geneList <- sort(geneList, decreasing = T)
        }
        
        fgseaRes <- fgsea(pathways = GS, 
                          stats = geneList,
                          minSize=15,
                          maxSize=500,
                          nperm=10000,
                          nproc=1)
        lead.edge <- fgseaRes[["leadingEdge"]]
        lead.edge <- lapply(lead.edge, function(x){ 
            x <- paste(x, collapse=",")
        })
        lead.edge <- unlist(lead.edge)
        
        fgseaRes[["lead.edge"]] <- lead.edge
        fgseaRes[["leadingEdge"]] <- NULL
        
        gseaList[[classes[j]]] <- setDF(fgseaRes)
    }
    
    list.es <- lapply(gseaList, "[", c("pathway", "ES"))
    list.es <- lapply(names(list.es), function(x) setNames(list.es[[x]], paste0(x, "_", names(list.es[[x]]))))
    es.combined <- Reduce(function(x, y) merge(x, y, all=T, by.x=grep("pathway", colnames(x)), by.y=grep("pathway", colnames(y))), list.es, accumulate=F)
    names(es.combined)[1] <- "pathway"
    
    #write.table(es.combined, file=paste0(name_out, "_Enrichment.txt"), sep="\t", row.names = F)
    
    list.pval <- lapply(gseaList, "[", c("pathway", "padj"))
    list.pval <- lapply(names(list.pval), function(x) setNames(list.pval[[x]], paste0(x, "_", names(list.pval[[x]]))))
    pval.combined <- Reduce(function(x, y) merge(x, y, all=T, by.x=grep("pathway", colnames(x)), by.y=grep("pathway", colnames(y))), list.pval, accumulate=F)
    names(pval.combined)[1] <- "pathway"
    
    #write.table(pval.combined, file=paste0(name_out, "_Enrichment_Pvalue.txt"), sep="\t", row.names = F)
    
    list.nes <- lapply(gseaList, "[", c("pathway", "NES"))
    list.nes <- lapply(names(list.nes), function(x) setNames(list.nes[[x]], paste0(x, "_", names(list.nes[[x]]))))
    nes.combined <- Reduce(function(x, y) merge(x, y, all=T, by.x=grep("pathway", colnames(x)), by.y=grep("pathway", colnames(y))), list.nes, accumulate=F)
    names(nes.combined)[1] <- "pathway"
    
    #write.table(nes.combined, file=paste0(name_out, "_Enrichment_NES.txt"), sep="\t", row.names = F)
    gsea.res <- list(es.combined, pval.combined, nes.combined)
    names(gsea.res) <- c("ES", "padj", "NES")
	message("Function doFastGSEA done.")
    return(gsea.res)
}    

suppressMessages(library('fgsea'))
suppressMessages(library('dplyr'))
suppressMessages(library('readr'))
suppressMessages(library('data.table'))

if (!interactive() && !exists('SOURCE')) {
	# Get and check arguments.
	suppressMessages(library(docopt))
	arg <- docopt(doc, version="0.0.1\n", strict=T)
	arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
	clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
	names(arg) <- clean(names(arg))

	# parameters
	#arg <- list()
	#arg$expression <- "tmp/cocozoes/GSE13052.tsv"
	#arg$sample_annotation <- "config/sample_annotation/GSE13052.tsv"
	#arg$annotation_cols <- c("ensembl_gene_id", "hgnc_symbol", "external_gene_name")
	#arg$symbols <- "external_gene_name"
	#arg$sample_name_col <- "Sample_geo_accession"
	#arg$class_col <- "Class"
	#arg$gmt <- "config/pathways/BTM.gmt"
	#arg$es_file <- "tmp/coco_es.tsv"
	#arg$nes_file <- "tmp/coco_nes.tsv"
	#arg$pval_file <- "tmp/coco_pval.tsv"


	exp.tibble <- read_tsv(arg$expression) %>%
		filter_(paste0("!is.na(", arg$symbols,")"))

	exp.df <- exp.tibble %>% 
		select_(.dots=paste0("-", arg$annotation_cols)) %>%
			as.data.frame
	rownames(exp.df) <- exp.tibble[[arg$symbols]]

	annot.df <- as.data.frame(read_tsv(arg$sample_annotation))
	rownames(annot.df) <- annot.df[[arg$sample_name_col]]

	gmt <- read.gmt(arg$gmt)
	geneset <- gmt[["genes"]]

	res <- doFastGSEA(exp.df, annot.df, geneset, class_col=arg$class_col)

	write_tsv(res[["ES"]], arg$es)
	write_tsv(res[["NES"]], arg$nes)
	write_tsv(res[["padj"]], arg$pval)
}
