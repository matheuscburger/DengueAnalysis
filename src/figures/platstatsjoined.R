#!/usr/bin/env Rscript

library("readr")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("ggrepel")
library("Vennerable")
library("gridExtra")
library("grid")


# Diretorio em que as figuras serao salvas
dir.create(file.path("figures", "platstatsjoined"))

# Le tabela de anotação das probes
annot <- read_tsv("config/reannotation/annotation_long.tsv")

# todas as plataformas
counts <- annot %>% filter(Hits == 1 & NumAnnot == 1) %>% 
	#filter(Platform %in% platforms) %>% 
    dplyr::select(Platform, Gene, Type, Group, Database) %>%
    unique() %>%
	count(Gene, Type, Group, Database) %>%
	mutate(Plats=n) %>% ungroup() %>%  
	dplyr::select(Plats, Group) %>% 
	count(Plats, Group) %>%
    group_by(Group) %>%
	arrange(desc(Plats))  %>%
	mutate(cumn=cumsum(n))

pl <- counts %>% filter(Group %in% c("coding", "lnoncoding")) %>% 
    ggplot(aes(x=factor(Plats), y=cumn, group=Group, color=Group)) +
        geom_point() +
	    geom_line(lwd=1.5) + 
		scale_color_hue(l=40) +
		theme_bw() +
        labs(x="Number of Platforms", y="Number of Genes", title="All platforms")

for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "platstatsjoined", paste0("all_plats.", ext)), plot=pl)
}
write_tsv(counts, file.path("figures", "platstatsjoined", "all_plats.counts.tsv"))

# minhas plataformas
platforms <- c("GPL570", "GPL2700")
counts <- annot %>% filter(Hits == 1 & NumAnnot == 1) %>% 
	filter(Platform %in% platforms) %>% 
    dplyr::select(Platform, Gene, Type, Group, Database) %>%
    unique() %>%
	count(Gene, Type, Group, Database) %>%
	mutate(Plats=n) %>% 
	ungroup() %>% 
	dplyr::select(Plats, Group) %>%
	count(Plats, Group) %>%
    group_by(Group) %>%
	arrange(desc(Plats)) %>%
	mutate(cumn=cumsum(n))

pl <- counts %>% filter(Group %in% c("coding", "lnoncoding")) %>% 
    ggplot(aes(x=factor(Plats), y=cumn, group=Group, color=Group)) +
        geom_point() +
		geom_text_repel(aes(label=cumn)) +
	    geom_line(lwd=1.5) + 
		scale_color_hue(l=40) +
		theme_bw() +
        labs(x="Number of Platforms", y="Number of Genes", title="GPL2700 and GPL570")

for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "platstatsjoined", paste0("my_plats.", ext)), plot=pl)
}
write_tsv(counts, file.path("figures", "platstatsjoined", "my_plats.counts.tsv"))

## Venn diagram

do_venn <- function(x){
	plat2genes <- lapply(platforms, function(p) dplyr::filter(x, Platform==p) %>% .[["Gene"]])
	names(plat2genes) <- platforms
	venn <- Venn(plat2genes)
	return(venn)
}

venn_list <- annot %>% filter(Hits == 1 & NumAnnot == 1) %>%  
	filter(Platform %in% platforms) %>% 
    dplyr::select(Platform, Gene, Group) %>% 
    unique() %>% group_by(Group) %>% 
	do(venn=do_venn(.))

fill_colors <- c("11"="#006694", "10"="#84D2F5", "01"="#84D2F5")

pdf(file.path("figures", "platstatsjoined", "venn_diagram.pdf"))
for(i in 1:nrow(venn_list)){
	x <- venn_list[i, ]
	plot_title <- x[["Group"]]
	venn <- x$venn[[1]]
	if(!is.null(venn) && sum(venn@IndicatorWeight[,'.Weight'] != 0) > 1){
		gp <- VennThemes(compute.Venn(venn))
		for(n in names(gp$SetText)) gp$SetText[[n]][["lwd"]] <- 0
		for(n in names(gp$Set)) gp$Set[[n]][["lwd"]] <- 0
		for(n in names(gp$Set)) gp$Set[[n]][["col"]] <- "white"
		for(n in names(gp$SetText)) gp$SetText[[n]][["col"]] <- "black"
		for(n in names(gp$Face)) gp$Face[[n]][["fill"]] <- fill_colors[n] 
		for(n in names(gp$FaceText)) gp$FaceText[[n]][["col"]] <- "white"
		p1 <- grid.grabExpr(plot(venn, gp=gp))
		grid.arrange(p1, top=plot_title)
	}
}
dev.off()
