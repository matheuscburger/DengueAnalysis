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

fill_colors2 <- c("11"="#006694", "10"="#84D2F5", "01"="#84D2F5")
fill_colors3 <- c("000"="#FFFFFF",
                  "100"="#84D2F5", "010"="#84D2F5", "001"="#84D2F5",
                  "110"="#00B1FF", "101"="#00B1FF", "100"="#00B1FF", "011"="#00B1FF", "010"="#00B1FF", "001"="#00B1FF",
                  "111"="#006694")
fill_colors4 <- c("0000"="#FFFFFF",
                  "1000"="#84D2F5", "0100"="#84D2F5", "0010"="#84D2F5", "0001"="#84D2F5",
                  "1100"="#00B1FF", "1010"="#00B1FF", "1001"="#00B1FF", "0110"="#00B1FF", "0101"="#00B1FF", "0011"="#00B1FF",
                  "0111"="#006694", "1110"="#006694", "1011"="#006694", "1101"="#006694", 
                  "1111"="#00374F")
pdf(file.path("figures", "venn", "venn.pdf"))
for(i in 1:nrow(deg_list)){
	x <- deg_list[["venn"]][[i]]
	plot_title <- deg_list[["title"]][[i]]
	if(!is.null(x$venn)){
		gp <- VennThemes(compute.Venn(x$venn))
		for(n in names(gp$SetText)) gp$SetText[[n]][["lwd"]] <- 0
		for(n in names(gp$Set)) gp$Set[[n]][["lwd"]] <- 0
		for(n in names(gp$Set)) gp$Set[[n]][["col"]] <- "white"
		for(n in names(gp$SetText)) gp$SetText[[n]][["col"]] <- "black"
		for(n in names(gp$FaceText)) gp$FaceText[[n]][["col"]] <- "white"

        if(ncol(x$venn@IndicatorWeight)-1 == 2){
            for(n in names(gp$Face)) gp$Face[[n]][["fill"]] <- fill_colors2[n] 
        }else if(ncol(x$venn@IndicatorWeight)-1 == 3){
            for(n in names(gp$Face)) gp$Face[[n]][["fill"]] <- fill_colors3[n] 
        }else if(ncol(x$venn@IndicatorWeight)-1 == 4){
            for(n in names(gp$Face)) gp$Face[[n]][["fill"]] <- fill_colors4[n] 
        }

		if(ncol(x$venn@IndicatorWeight)-1 > 3){
			p1 <- grid.grabExpr(plot(x$venn, gp=gp, type="ellipses"))
		} else {
			p1 <- grid.grabExpr(plot(x$venn, gp=gp))
		}
		grid.arrange(p1, top=plot_title)
	}
}
dev.off()


