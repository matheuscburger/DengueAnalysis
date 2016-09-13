
library("readr")
library("dplyr")
library("ggplot2")
library("ggthemes")


# Diretorio em que as figuras serao salvas
dir.create(file.path("figures", "platstatsjoined"))

# Le tabela de anotação das probes
annot <- read_tsv("config/reannotation/annotation_long.tsv")

# todas as plataformas
counts <- annot %>% filter(Hits == 1 & NumAnnot == 1) %>% 
	#filter(Platform %in% platforms) %>% 
    select(Platform, Gene, Type, Group, Database) %>%
    unique() %>%
	count(Gene, Type, Group, Database) %>%
	mutate(Plats=n) %>% 
	count(Plats, Group) %>%
    group_by(Group) %>%
	arrange(desc(Plats)) %>%
	mutate(cumn=cumsum(n))

pl <- counts %>% filter(Group %in% c("coding", "lnoncoding")) %>% 
    ggplot(aes(x=factor(Plats), y=cumn, group=Group, color=Group)) +
        geom_point() +
	    geom_line(lwd=1.5) + 
		scale_color_hue(l=40) +
		theme_bw() +
        labs(x="Number of Platforms", y="Number of Genes")

for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "platstatsjoined", paste0("all_plats.", ext)), plot=pl)
}

# minhas plataformas
platforms <- c("GPL570", "GPL2700")
counts <- annot %>% filter(Hits == 1 & NumAnnot == 1) %>% 
	filter(Platform %in% platforms) %>% 
    select(Platform, Gene, Type, Group, Database) %>%
    unique() %>%
	count(Gene, Type, Group, Database) %>%
	mutate(Plats=n) %>% 
	count(Plats, Group) %>%
    group_by(Group) %>%
	arrange(desc(Plats)) %>%
	mutate(cumn=cumsum(n))

pl <- counts %>% filter(Group %in% c("coding", "lnoncoding")) %>% 
    ggplot(aes(x=factor(Plats), y=cumn, group=Group, color=Group)) +
        geom_point() +
	    geom_line(lwd=1.5) + 
		scale_color_hue(l=40) +
		theme_bw() +
        labs(x="Number of Platforms", y="Number of Genes")

for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "platstatsjoined", paste0("my_plats.", ext)), plot=pl)
}
