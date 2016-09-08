# Fazer um gráfico de pizz com as classes das amostras em cada estudo (coluna ExtendedClass e Class)

options(stringsAsFactors=FALSE)

library("readr")
library("tibble")
library("tidyr")
library("dplyr")
library("ggplot2")

# Diretorio em que as figuras serao salvas
dir.create(file.path("figures", "StudyInformation"))

# Diretorio com anotcao das amostras
sa_dir <- "config/sample_annotation/"

# Estudos para olhar
studies <- c("GSE13052", "GSE28405", "GSE43777", "GSE51808")

# Le anotacao das amostras, pega colunas de interesse e combina estudos
all_df <- data.frame()
for( s in studies ){
	s_df <- read_tsv(file.path(sa_dir, paste0(s, ".tsv")))
	all_df <- s_df %>% select(Sample_series_id, Sample_platform_id, Sample_geo_accession, ExtendedClass, Class) %>% rbind(all_df)
}

# Remove plataforma GPL201
all_df <- all_df %>% filter(Sample_platform_id != "GPL201")

# Faz pie chart para ExtendedClass

extendedClass_levels <- c("DSS", "DHF", "DF", "Dengue",
						  "Convalescent", "Control", "Healthy")
# obtem informacao resumida 
sum_info_extclass <- all_df %>% count(Sample_series_id, ExtendedClass) %>% 
	group_by(Sample_series_id) %>% 
	mutate(prop=n/sum(n), 
		   ExtendedClass = factor(ExtendedClass, levels=extendedClass_levels)) %>%
	mutate(pos=cumsum(prop) - 0.5*prop)

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
# faz grafico
extendedClass_colors <- c( # mais quente 
						  "DSS" = "#9F2110",
						  "DHF" = "#D76136",
						  "DF" = "#E5D43E",
						  "Dengue" = "#E59C00",
						  # mais frio
						  "Convalescent" = "#00C2D2",
						  "Control" = "#304D78",
						  "Healthy" = "#164D9F")
pl <- ggplot(sum_info_extclass, aes(x=factor(1), y=n, fill=ExtendedClass)) +
	geom_bar(stat="identity", position="fill", width=1) +
	geom_text(aes(label=n, y=pos), color="white", size=5, fontface="bold") +
	coord_polar(theta="y") +
	facet_grid(. ~ Sample_series_id) +
	scale_fill_manual(values=extendedClass_colors) +
	labs(fill="Extended Class", title="Studies") + theme_simple
for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "StudyInformation", paste0("extended_class.", ext)), plot=pl, width=10, height=3)
}


# Faz pie chart para Class

class_levels <- c("Control", "Dengue")
# obtem informacao resumida 
sum_info_class <- all_df %>% count(Sample_series_id, Class) %>% 
	group_by(Sample_series_id) %>% 
	mutate(prop=n/sum(n), 
		   Class = factor(Class, levels=class_levels)) %>%
	mutate(pos=cumsum(prop) - 0.5*prop)

# faz grafico
pl <- ggplot(sum_info_class, aes(x=factor(1), y=n, fill=Class)) +
	geom_bar(stat="identity", position="fill", width=1) +
	geom_text(aes(label=n, y=pos), color="white", size=5, fontface="bold") +
	coord_polar(theta="y") +
	facet_grid(. ~ Sample_series_id) +
	scale_fill_manual(values=extendedClass_colors) +
	labs(fill="Class", title="Studies") + theme_simple
for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "StudyInformation", paste0("class.", ext)), plot=pl, width=10, height=3)
}
