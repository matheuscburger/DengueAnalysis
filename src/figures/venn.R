#!/usr/bin/env Rscript

library("dplyr")
library("readr")
library("tidyr")
library("Vennerable")
library("gridExtra")
library("grid")

# setwd("../..")

# Diretorio em que as figuras serao salvas
dir.create(file.path("figures", "venn"))

# Carregar lncRNAs DEGs
degs <- read_tsv("results/joined_DEG.tsv")

do_venn <- function(x) {
	cols <- grep("GSE", names(x), value=T) 
	groups <- lapply(cols, function(y) x[["Symbol"]][grep("(Down)|(Up)", x[[y]])] )
	names(groups) <- cols
	rm_groups <- -which(sapply(groups, length) == 0)
	if(length(rm_groups) > 0){
		notempty <- groups[-which(sapply(groups, length) == 0)]
	} else {
		notempty <- groups
	}
	if(length(notempty) >= 2){
		venn <- Venn(notempty)
	} else {
		venn <- NULL
	}
	return(list(venn=venn, groups=groups))
}
deg_list <- degs %>% group_by(title) %>% do(venn=do_venn(.)) 

pdf(file.path("figures", "venn", "venn.pdf"))
for(i in 1:nrow(deg_list)){
	x <- deg_list[["venn"]][[i]]
	plot_title <- deg_list[["title"]][[i]]
	if(!is.null(x$venn)){
		gp <- VennThemes(compute.Venn(x$venn))
		for(n in names(gp$SetText)) gp$SetText[[n]][["lwd"]] <- 0
		for(n in names(gp$Set)) gp$Set[[n]][["lwd"]] <- 0
		if(ncol(x$venn@IndicatorWeight)-1 > 3){
			p1 <- grid.grabExpr(plot(x$venn, gp=gp, type="ellipses"))
		} else {
			p1 <- grid.grabExpr(plot(x$venn, gp=gp))
		}
		grid.arrange(p1, top=plot_title)
	}
}
dev.off()


