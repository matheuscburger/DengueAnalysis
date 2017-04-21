#!/usr/bin/env Rscript

library("readr")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")


# Diretorio em que as figuras serao salvas
dir.create(file.path("figures", "platstats"))

# Le tabela de anotação das probes
annot <- read_tsv("config/reannotation/annotation_long.tsv")
##################################################################################
# 1) Grafico de pizzas mostrando numero de probes por database em cada plataforma

# Pega informacao de databases por plataforma
db_levels <- annot %>% dplyr::select(Database) %>% table() %>% sort(decreasing=TRUE) %>% names
db_count <- annot %>% count(Platform, Database) %>% 
	group_by(Platform) %>% 
	mutate(prop=n/sum(n),
		   Database=factor(Database, level=db_levels)) %>%
    arrange(Platform, desc(Database)) %>%
	mutate(pos=cumsum(prop) - 0.5*prop)

# simplifica tema para grafico
theme_simple <- theme_minimal() +
	theme(axis.line=element_blank(),
		  axis.text.x=element_blank(),
          axis.text.y=element_blank(),
		  axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
 		  panel.background=element_blank(),
		  panel.border=element_blank(),
		  panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
		  plot.background=element_blank() 
		  )

# Faz grafico


database_colors <- c("gencode"="#E59C00",
					 "noncode"="#00A85E",
					 "lncipedia"="#8712B5",
					 "mitranscriptome"="#D9CA00",
					 "unannotated"="#304D78")
pl <- ggplot(db_count, aes(x=factor(1), y=n, fill=Database)) +
	geom_bar(stat="identity", position="fill", width=1) +
	geom_text(aes(label=n, y=pos), color="white", size=3, fontface="bold") +
	coord_polar(theta="y") +
	facet_wrap(~Platform) +
	scale_fill_manual(values=database_colors) +
	labs(fill="Databases", title="Number of Probes by Databases") + 
	theme_simple + theme(text = element_text(size=20))
for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "platstats", paste0("databases.", ext)), plot=pl)
}
write_tsv(db_count, file.path("figures", "platstats", paste0("databases.db_count.tsv")))


##################################################################################
# 2) Grafico de barras mostrando números de probes em cada grupo

change_colors <- function(incolors){
	myhsv <- rgb2hsv(col2rgb(incolors))
	myhsv["s", ] <- sapply(myhsv["s", ] + 0.8, min, 1)
	myhsv["v", ] <- sapply(myhsv["v", ] - 0.4, max, 0)
	res <- apply(myhsv, 2, function(x) hsv(h=x[1], s=x[2], v=x[3]))
	return(res)
}
# pegar os caras que sao multi hits e multi annotations
group_levels <- annot %>% filter(NumAnnot == 1 & Hits ==1) %>% dplyr::select(Group) %>% table() %>%
   	sort(decreasing=TRUE) %>% names %>%
	c("MultipleHits", "MultipleAnnotations")
group_colors <- change_colors(brewer.pal(length(group_levels), "Set1"))
names(group_colors) <- group_levels


group_count <- annot %>% 
    mutate(Group=replace(Group, NumAnnot > 1, "MultipleAnnotations")) %>%
	mutate(Group=replace(Group, Hits > 1, "MultipleHits")) %>% 
    count(Platform, Group) %>%	group_by(Platform) %>% 
	mutate(Group=factor(Group, level=group_levels)) %>%
	mutate(pos=0.5*n)

max_n <- max(group_count[, "n"])


pdf(file=file.path("figures", "platstats", paste0("groups.pdf")), width=10)
for(plat in unique(group_count[["Platform"]])){
	pl <- group_count %>% filter(Platform == plat) %>%
		ggplot(aes(x=Group, y=n)) +
		geom_bar(stat="identity", aes(fill=Group)) +
	    geom_text(aes(label=n, y=pos), color="white", size=5, fontface="bold") +
		ylim(0, max_n) +
		scale_fill_manual(values=group_colors) +
		labs(y="Count", x="Group", title=plat) +
		theme_bw() +
		theme(axis.text.x = element_text(angle=45, hjust=1)) +
	    theme(text = element_text(size=20))
	print(pl)
}
dev.off()
write_tsv(group_count, file.path("figures", "platstats", "groups.group_count.tsv"))

##################################################################################
# 3) Grafico de barras mostrando números de probes em cada tipo
type_levels <- annot %>% filter(NumAnnot == 1 & Hits == 1 & Group == "lnoncoding") %>%
   	dplyr::select(Type) %>% table() %>%
   	sort(decreasing=TRUE) %>% names %>%
	c("MultipleHits", "MultipleAnnotations")
type_colors <- change_colors(brewer.pal(length(type_levels), "Set3"))
names(type_colors) <- type_levels
type_count <- annot %>% filter(Group == "lnoncoding") %>% 
    mutate(Type=replace(Type, NumAnnot > 1, "MultipleAnnotations")) %>%
	mutate(Type=replace(Type, Hits > 1, "MultipleHits")) %>% 
    count(Platform, Type) %>% group_by(Platform) %>% 
	mutate(Type=factor(Type, level=type_levels)) %>%
	mutate(pos=0.5*n)

max_n <- max(type_count[, "n"])

pdf(file=file.path("figures", "platstats", paste0("types.pdf")), width=10)
for(plat in unique(group_count[["Platform"]])){
	pl <- type_count %>% filter(Platform == plat) %>%
		ggplot(aes(x=Type, y=n)) +
		geom_bar(stat="identity", aes(fill=Type)) +
		ylim(0, max_n) +
		scale_fill_manual(values=type_colors) +
	    geom_text(aes(label=n, y=pos), color="white", size=5, fontface="bold") +
		labs(y="Count", x="Type", title=plat) +
		theme_bw() +
		theme(axis.text.x = element_text(angle=45, hjust=1)) +
	    theme(text = element_text(size=20))
	print(pl)
}
dev.off()

write_tsv(type_count, file.path("figures", "platstats", "types.type_count.tsv"))
