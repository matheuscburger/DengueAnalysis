
library("readr")
library("dplyr")
library("ggplot2")
library("ggthemes")


# Diretorio em que as figuras serao salvas
dir.create(file.path("figures", "platstats"))

# Le tabela de anotação das probes
annot <- read_tsv("config/reannotation/annotation_long.tsv")
##################################################################################
# 1) Grafico de pizzas mostrando numero de probes por database em cada plataforma

# Pega informacao de databases por plataforma
db_levels <- annot %>% select(Database) %>% table() %>% sort(decreasing=TRUE) %>% names
db_count <- annot %>% count(Platform, Database) %>% 
	group_by(Platform) %>% 
	mutate(prop=n/sum(n),
		   Database=factor(Database, level=db_levels)) %>%
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
	theme_simple
for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "platstats", paste0("databases.", ext)), plot=pl)
}


##################################################################################
# 2) Grafico de barras mostrando números de probes em cada grupo

# pegar os caras que sao multi hits e multi annotations
group_levels <- annot %>% select(Group) %>% table() %>%
   	sort(decreasing=TRUE) %>% names %>%
	c("MultipleHits", "MultipleAnnotations")
group_count <- annot %>% 
    mutate(Group=replace(Group, NumAnnot > 1, "MultipleAnnotations")) %>%
	mutate(Group=replace(Group, Hits > 1, "MultipleHits")) %>% 
    count(Platform, Group) %>%	group_by(Platform) %>% 
	mutate(prop=n/sum(n),
		   Group=factor(Group, level=group_levels)) %>%
	mutate(pos=cumsum(prop) - 0.5*prop)

max_n <- max(group_count[, "n"])


pdf(file=file.path("figures", "platstats", paste0("groups.pdf")))
for(plat in unique(group_count[["Platform"]])){
	pl <- group_count %>% filter(Platform == plat) %>%
		ggplot(aes(x=Group, y=n)) +
		geom_bar(stat="identity", aes(fill=Group)) +
		ylim(0, max_n) +
		scale_fill_hue(l=40) +
		labs(y="Count", x="Group", title=plat) +
		theme_bw() +
		theme(axis.text.x = element_text(angle=45, hjust=1))
	print(pl)
}
dev.off()

##################################################################################
# 3) Grafico de barras mostrando números de probes em cada tipo
type_levels <- annot %>% select(Type) %>% table() %>%
   	sort(decreasing=TRUE) %>% names %>%
	c("MultipleHits", "MultipleAnnotations")
type_count <- annot %>% filter(Group == "lnoncoding") %>% 
    mutate(Type=replace(Type, NumAnnot > 1, "MultipleAnnotations")) %>%
	mutate(Type=replace(Type, Hits > 1, "MultipleHits")) %>% 
    count(Platform, Type) %>% group_by(Platform) %>% 
	mutate(prop=n/sum(n),
		   Type=factor(Type, level=type_levels)) %>%
	mutate(pos=cumsum(prop) - 0.5*prop)

max_n <- max(type_count[, "n"])

pdf(file=file.path("figures", "platstats", paste0("types.pdf")))
for(plat in unique(group_count[["Platform"]])){
	pl <- type_count %>% filter(Platform == plat) %>%
		ggplot(aes(x=Type, y=n)) +
		geom_bar(stat="identity", aes(fill=Type)) +
		#ylim(0, max_n) +
		scale_fill_hue(l=40) +
		labs(y="Count", x="Type", title=plat) +
		theme_bw() +
		theme(axis.text.x = element_text(angle=45, hjust=1))
	print(pl)
}
dev.off()
