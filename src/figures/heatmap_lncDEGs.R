#!/usr/bin/env Rscript

# Carregar bibliotecas
library("readr")
library("dplyr")
library("tidyr")
library("stringr")
library("Heatplus")
library("pvclust")
library("dendextend")
library("reshape2")


# Diretorio em que as figuras serao salvas
dir.create(file.path("figures", "heatmap_lncDEGs"))

# Carregar lncRNAs DEGs
degs <- read_tsv("results/joined_DEG.tsv")
lnc.tbl <- degs %>% 
	filter(title %in% c("Dengue vs Control", "EarlyDengue vs Control", "LateDengue vs Control")) %>%
	filter(Ratio >= 0.6) %>%
	filter(!HasDiscordantDEGs) %>%
	filter(Group == "lnoncoding")
lncs <- unique(lnc.tbl[['Symbol']])

comps <- lnc.tbl %>% dplyr::select(Symbol, title) %>% 
	spread(key=title, value=title)
inds <- which(!is.na(comps[, 2:ncol(comps)]), arr.ind=T)
inds[,2] <- inds[,2] + 1
comps[inds] <- TRUE
comps[which(is.na(comps), arr.ind=T)] <- FALSE

# Carregar arquivos de expressão
studies <- c("GSE13052", "GSE28405", "GSE43777", "GSE51808")
exp.data.dir <- "data/processed/filtered/"
files <- paste0(exp.data.dir, studies, ".tsv")
exps <- lapply(files, read_tsv)
names(exps) <- studies

# pega somente lncs em exps
exps.lncs <- lapply(exps, function(x) x %>% filter(Symbol %in% lncs))

# Carregar arquivos de anotacao das amostras
sannot.data.dir <- "config/sample_annotation/"
files <- paste0(sannot.data.dir, studies, ".tsv")
sannot <- lapply(files, read_tsv)
names(sannot) <- studies

# Fazer Heatmap
rb <- colorRampPalette(c(rgb(0,0,1), rgb(1, 1, 1), rgb(1,0,0)), alpha = TRUE)(100)
my.hclust <- function(x) hclust(x, method="ward.D2")
dist.euclidean <- function(x) dist(x, method="euclidean")

for( i in 1:4 ){
	pdf(paste0("figures/heatmap_lncDEGs/", studies[i], ".pdf"))
	curr.exps <- exps.lncs[[i]]
	curr.sannot <- sannot[[i]]
	curr.samples <- curr.sannot[["Sample_geo_accession"]]
	curr.comps <- full_join(curr.exps[, c("ProbeName", "Symbol")], comps) %>%
		slice(match(curr.exps[["ProbeName"]], ProbeName)) %>% 
		dplyr::select(-ProbeName, -Symbol) %>% apply(2, as.logical)


	scaled <- t(apply(curr.exps[, curr.samples], 1, scale))
	colnames(scaled) <- curr.samples

	#hcrow <- hclust(dist(scaled, method="euclidean"), method="ward.D2")
	hccol <- hclust(dist(t(scaled), method="euclidean"), method="ward.D2")
	cutheight <- sort(get_nodes_attr(as.dendrogram(hccol), "height"), decreasing=T)[2] - 1

	#heat <- annHeatmap(as.matrix(curr.exps[, curr.samples]), 
	heat <- annHeatmap2(scaled, 
					   col=rb, breaks=seq(-3, 3, length=101),
					   annotation=list(Col=list(data=curr.sannot[, c("Class", "Stage")])),
									   #Row=list(data=curr.comps)),
					   labels=list(Row=list(labels=curr.exps[["Symbol"]], nrow=5, cex=.6),
								   Col=list(nrow=5)),
					   scale="none",
					   cluster=list(cuth=cutheight),
					   dendrogram=list(clustfun=my.hclust, distfun=dist.euclidean),
					   legend=T)
	plot(heat)
	dev.off()
}
