#!/usr/bin/env Rscript
# exp out outdir mode pV [temp] [pA]

"
Usage:
	genefilter.R --expr-file=<file> --output=<file> --outdir=<dir> --method=<method> --pvalue=<value> [--template=<file> --pvanova=<value>]

Options:
	-h --help	show this help message
	--version	show program version
	--expr-file=<file>	expression file
	--output=<file>	output file
	--outdir=<dir>	output directory
	--method=<method>	filtering method
	--pvalue=<value>	pvalue for variance
	--template=<file>	sample annotation file for ANOVA
	--pvanova=<value>	pvalue for ANOVA
" -> doc

library("docopt")
library("pracma")
library("data.table")
options(stringsAsFactors=FALSE)
args <- docopt(doc, version="0.0.1\n", strict=T)
args <- args[!sapply(args, is.null)]
clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
names(args) <- clean(names(args))
print(args)
#args <- commandArgs(TRUE)

doVar <- function(D){
	pV <<- as.numeric(args[["pvalue"]])
	V <- apply(D, 1, var)
	length(V)

	E <- mean(V)
	E2 <- mean(V**2)
	
	ah <- E**2/(E2 - E**2) + 2
	bh <- (ah-1)*(ah-2)*(E2 - E**2)/E
	
	Phist <- function(x){
		IG <- gammainc(bh/x, ah)["uppinc"]
		G <- gamma(ah)
		return(1 - IG/G)
	}
	
	P <- sapply(V, Phist)
	names(P) <- gsub(".uppinc", "", names(P))
	length(P)
	return(P)
}


#doANOVA <- function(D){
#	if(length(args) != 7) stop("Hey! You forgot something!")
#	pA <<- as.numeric(args[["pvanova"]])
#	temp <- read.delim(args[["template"]], header=TRUE, row.names=1)
#	Pval <- apply(D[, rownames(temp)], 1, function(x){
#		aov.res <- aov(x ~ temp$Class)
#		return(summary(aov.res)[[1]][["Pr(>F)"]][1])
#	})
#	if(args[["method"]] == "ANOVA") Pval <- p.adjust(Pval, method="fdr")
#	return(Pval)
#}
#
#r1 <- function(x){
#	m = (pA - 1.)/pV
#	return(pA + m*(x - pV))
#}
#r2 <- function(x){
#	m = -pA/(1. - pV)
#	return(pA + m*(x - pV))
#}

main <- function(){
	# all praise to the lord data.table
	D0 <- fread(args[["expr_file"]], header=T, data.table=F)
    rown <- D0[,1]
    D0 <- D0[,-1]
	rownames(D0) <- rown	
	PV <- 0
    PAoV <- 0
    	
	#D0 <- read.table(args[1], header=TRUE, row.names=1, sep="\t")
	initialdim <- dim(D0)
	print(paste0('TAMANHO INICIAL: ', initialdim))
	E <- apply(D0, 1, mean)
	D0 <- D0[order(E, decreasing=TRUE), ]
	D0 <- D0[1:(.8*nrow(D0)),]

	if(args[["method"]] == "V"){
		PV <- doVar(D0)
		print(paste0('length PV: ', length(PV)))
		DEGs <- names(PV)[PV < pV]
		print(paste0('length DEGs: ', length(DEGs)))
	}
#	else if(args[["method"]] == "ANOVA"){
#		PAoV <- doANOVA(D0)
#		DEGs <- names(PAoV)[PAoV < pA]
#	}
#	else if(args[["method"]] == "AV"){
#		PV <- doVar(D0)
#		PAoV <- doANOVA(D0)
#	
#		PP <- rbind(PV, PAoV)
#	
#		DEGs <- which(PP[2, ] < r1(PP[1, ]) | PP[2, ] < r2(PP[1, ]))
#	}
#	else{
#		stop("WHAT ARE YOU DOING?? I've never heard of this method!")
#	}
	

	D0 <- D0[DEGs, ]
	print(paste0('dim DEGs: ', dim(D0)))
	
	D0 <- cbind(Symbol=rownames(D0), D0)
    filterdim <- dim(D0)
	
	outdir <- paste0(args[["outdir"]],"/")
	write.table(D0, paste0(outdir, args[["output"]]),
					sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}


main()
