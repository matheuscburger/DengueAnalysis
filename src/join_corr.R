#!/usr/bin/env Rscript
SOURCE <- TRUE
source("src/microarrayAnalysis/join_stats.R")
source("src/microarrayAnalysis/collapse.R")


filesGPL570 <- c("results/correlation/byprobe/GSE43777.tsv", "results/correlation/byprobe/GSE51808.tsv")
filesGPL2700 <- c("results/correlation/byprobe/GSE13052.tsv", "results/correlation/byprobe/GSE28405.tsv")

library("stringr")
library("data.table")

maintain.cols <- c("correlation", "p.value", "adj.p.value")
key.cols <- c("mrn_id", "lnc_id", "int_type")
probe.col <- c("ProbeName_mrn", "ProbeName_lnc")

joinedGPL570 <- as.data.frame(join.stats(filesGPL570, maintain.cols=maintain.cols, 
				  key.cols=c(key.cols, probe.col),
				  probe.col=probe.col))
write.table(joinedGPL570, "results/correlation/GPL570.tsv", sep="\t", quote=FALSE, row.names=FALSE)

joinedGPL2700 <- as.data.frame(join.stats(filesGPL2700, maintain.cols=maintain.cols, 
				  key.cols=c(key.cols, probe.col),
				  probe.col=probe.col))
write.table(joinedGPL2700, "results/correlation/GPL2700.tsv", sep="\t", quote=FALSE, row.names=FALSE)

maxcor <- function(x){
	x[order(abs(x[, "correlation"]), decreasing=TRUE)[1], ]
}
collapse_cor <- function(input, output,  by_col){
	in.df <- read.delim(input)
	exp.df <- collapse(in.df, method=maxcor, by_col=by_col, drop=FALSE)
	write.table(exp.df, output, quote=FALSE, sep="\t", row.names=FALSE)
}


by_genes_files <- c()
for(f in c(filesGPL570, filesGPL2700)){
	base <- basename(f)
	res_fname <- file.path("results/correlation/bygene", base)
	by_genes_files <- c(by_genes_files, res_fname)
	collapse_cor(f, res_fname, by_col=key.cols)
}

joined <- as.data.frame(join.stats(by_genes_files, maintain.cols=c(maintain.cols, probe.col), 
				  key.cols=key.cols))


joined$min.cor <- apply(joined[, grep("^correlation", colnames(joined), value=T)], 1, function(x){
	if(any(is.na(x))){
		return(NA)
	} else {
		return(x[which.min(abs(x))])
	}
})

joined$max.p <- apply(joined[, grep("^p.value", colnames(joined), value=T)], 1, max)
joined$max.adj.p <- apply(joined[, grep("^adj.p.value", colnames(joined), value=T)], 1, max)

write.table(joined, "results/correlation/all.tsv", sep="\t", quote=FALSE, row.names=FALSE)
