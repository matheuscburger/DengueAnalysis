#!/usr/bin/env Rscript

library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggthemes")
library("pbapply")

# Diretorio em que as figuras serao salvas
dir.create(file.path("figures", "DEG_cutoff"))

# le anotacao
reannotation <- read_tsv("config/reannotation/annotation_long.tsv") %>% 
	dplyr::select(Database, Gene, Type, Group) %>% unique

# le estatisticas
gpl570 <- read_tsv("results/GPL570_joined_DEG.tsv", guess_max=100000)
gpl2700 <- read_tsv("results/GPL2700_joined_DEG.tsv", guess_max=10000000)

# conta o numero de estudos por teste
countStudies <- function(x){
	# pega Log2FC
	log2fccols <- grep("Log2FC", colnames(x), value=T)
	# quantas colunas de log2FC possuem NA ao inves de valores
	# para descobrir quais sao os estudos que possuem uma comparacao
	countl2fc <- apply(apply(x[, log2fccols], 2, is.na), 1, function(x) sum(!x))
	# pega o maximo para uma determinada comparacao
	nstudies <- sapply(split(countl2fc, x[, "title"]), max)
	x[, "Log2FC_Count"] <- countl2fc
	x[, "NStudies"] <- nstudies[x[["title"]]]
	x[, "HasNAs"] <- x[, "NStudies"] != x[, "Log2FC_Count"]
	return(x)
}

gpl570 <- countStudies(gpl570)
gpl2700 <- countStudies(gpl2700)

# verifica se fold-change eh discordante
is_discordant <- function(x){
	l <- length(unique(as.numeric(na.exclude(sign(x)))))
	if(l > 1){
		return(T)
	} else {
		return(F)
	}
}

# pega FCs
# GPL2700
message("O dataset GPL2700 possui ", nrow(gpl2700), " linhas")
sub_gpl2700 <- gpl2700 %>% filter(!HasNAs) %>% # remove linhas que possuem NAs devido a filtro 
	# Seleciona colunas de interesse
	dplyr::select(title, ProbeName, Symbol, starts_with("Log2FC"), starts_with("LIMMA.rawp"), starts_with("LIMMA.adjp")) %>% 
	left_join(reannotation, by=c("Symbol" = "Gene")) # Adiciona informacoes da reanotacao
message("O dataset GPL2700 possui ", nrow(sub_gpl2700), " linhas apos remocao de NAs")
# verifica fold-changes discordantes
sub_gpl2700[, "Discordant"] <- dplyr::select(sub_gpl2700, starts_with("Log2FC.GSE")) %>% apply(1, is_discordant)
# verifica menor fold-change
sub_gpl2700[, "min_fc"] <- dplyr::select(sub_gpl2700, starts_with("Log2FC.GSE")) %>% 
	apply(1, function(x) {min(abs(x), na.rm=T)})
# verifica maior p-valor
sub_gpl2700[, "max_p"] <- dplyr::select(sub_gpl2700, starts_with("LIMMA.rawp")) %>% 
	apply(1, function(x) {max(abs(x), na.rm=T)})
# verifica maior p-valor ajustado
sub_gpl2700[, "max_adjp"] <- dplyr::select(sub_gpl2700, starts_with("LIMMA.adjp")) %>% 
	apply(1, function(x) {max(abs(x), na.rm=T)})

# GPL570
message("O dataset GPL570 possui ", nrow(gpl570), " linhas")
sub_gpl570 <- gpl570 %>% filter(!HasNAs) %>%
	dplyr::select(title, ProbeName, Symbol, starts_with("Log2FC"), starts_with("LIMMA.rawp"), starts_with("LIMMA.adjp")) %>% 
	left_join(reannotation, by=c("Symbol" = "Gene"))
message("O dataset GPL570 possui ", nrow(sub_gpl570), " linhas apos remocao de NAs")
sub_gpl570[, "Discordant"] <- dplyr::select(sub_gpl570, starts_with("Log2FC.GSE")) %>% apply(1, is_discordant)
sub_gpl570[, "min_fc"] <- dplyr::select(sub_gpl570, starts_with("Log2FC.GSE")) %>% 
	apply(1, function(x) {min(abs(x), na.rm=T)})
sub_gpl570[, "max_p"] <- dplyr::select(sub_gpl570, starts_with("LIMMA.rawp")) %>% 
	apply(1, function(x) {max(abs(x), na.rm=T)})
sub_gpl570[, "max_adjp"] <- dplyr::select(sub_gpl570, starts_with("LIMMA.adjp")) %>% 
	apply(1, function(x) {max(abs(x), na.rm=T)})

# pega range de FCs
all_fc <- c(sub_gpl570 %>% dplyr::select(starts_with("Log2FC")) %>% unlist %>% as.numeric, 
			sub_gpl2700 %>% dplyr::select(starts_with("Log2FC")) %>% unlist %>% as.numeric)
max_fc <- max(abs(all_fc), na.rm=T)

# funcao para pegar tamanho da interseccao
get_intersection_size_fc <- function(x, comp){
	co <- x["cutoff"]
	symbols_gpl570 <- sub_gpl570 %>% filter(min_fc >= co, title == comp, Group == "lnoncoding", !Discordant) %>% .[["Symbol"]] 
	symbols_gpl2700 <- sub_gpl2700 %>% filter(min_fc >= co, title == comp, Group == "lnoncoding", !Discordant) %>% .[["Symbol"]] 
	symbols_alls <- intersect(unique(symbols_gpl570), unique(symbols_gpl2700))
	return(length(symbols_alls))
}

# calcula estatisticas 
cutoff.df <- data.frame(cutoff=seq(0, max_fc, length.out=1000))
cutoff.df[, "EarlyDengue"] <- pbapply(cutoff.df, 1, get_intersection_size_fc, comp="EarlyDengue vs Control")
cutoff.df[, "LateDengue"] <- pbapply(cutoff.df, 1, get_intersection_size_fc, comp="LateDengue vs Control")
cutoff.df[, "Dengue"] <- pbapply(cutoff.df, 1, get_intersection_size_fc, comp="Dengue vs Control")

cutoff.melted <- cutoff.df %>% gather("Comparison", "DEGs", c(EarlyDengue, LateDengue, Dengue)) %>% 
	filter(DEGs != 0)

(pl_fc <- ggplot(cutoff.melted, aes(x=cutoff, y=DEGs, color=Comparison)) +
	geom_line(size=1) +
	geom_vline(xintercept=log2(1.5), linetype = "longdash", color="gray") +
	geom_text(label="FC = 1.5", x=log2(1.5), y=Inf, vjust=4, hjust=-.1, color="darkgray") +
	geom_vline(xintercept=log2(1.25), linetype = "longdash", color="gray") +
	geom_text(label="FC = 1.25", x=log2(1.25), y=Inf, vjust=2, hjust=-.1, color="darkgray") +
	labs(x="Log2(Fold-Change)", title="Number of Lnoncoding vs Cutoff Used") +
	theme_bw())
for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "DEG_cutoff", paste0("fold_change.", ext)), plot=pl_fc)
}

# p-valores

# funcao para pegar tamanho da interseccao
get_intersection_size_pv <- function(x, comp){
	co <- x["cutoff"]
	symbols_gpl570 <- sub_gpl570 %>% filter(max_p <= co, title == comp, Group == "lnoncoding", !Discordant) %>% .[["Symbol"]] 
	symbols_gpl2700 <- sub_gpl2700 %>% filter(max_p <= co, title == comp, Group == "lnoncoding", !Discordant) %>% .[["Symbol"]] 
	symbols_alls <- intersect(unique(symbols_gpl570), unique(symbols_gpl2700))
	return(length(symbols_alls))
}

# calcula estatisticas 
cutoff.df <- data.frame(cutoff=seq(0, 0.25, length.out=1000))
cutoff.df[, "EarlyDengue"] <- pbapply(cutoff.df, 1, get_intersection_size_pv, comp="EarlyDengue vs Control")
cutoff.df[, "LateDengue"] <- pbapply(cutoff.df, 1, get_intersection_size_pv, comp="LateDengue vs Control")
cutoff.df[, "Dengue"] <- pbapply(cutoff.df, 1, get_intersection_size_pv, comp="Dengue vs Control")

cutoff.melted <- cutoff.df %>% gather("Comparison", "DEGs", c(EarlyDengue, LateDengue, Dengue)) %>% 
	filter(DEGs != 0)

(pl_pv <- ggplot(cutoff.melted, aes(x=cutoff, y=DEGs, color=Comparison)) +
	geom_line(size=1) +
	geom_vline(xintercept=0.05, linetype = "longdash", color="gray") +
	geom_vline(xintercept=0.01, linetype = "longdash", color="gray") +
	labs(x="P-value", title="Number of Lnoncoding vs Cutoff Used") +
	theme_bw())
for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "DEG_cutoff", paste0("pvalue.", ext)), plot=pl_pv)
}


# p-valores ajustados

# funcao para pegar tamanho da interseccao
get_intersection_size_adjpv <- function(x, comp){
	co <- x["cutoff"]
	symbols_gpl570 <- sub_gpl570 %>% filter(max_adjp <= co, title == comp, Group == "lnoncoding", !Discordant) %>% .[["Symbol"]] 
	symbols_gpl2700 <- sub_gpl2700 %>% filter(max_adjp <= co, title == comp, Group == "lnoncoding", !Discordant) %>% .[["Symbol"]] 
	symbols_alls <- intersect(unique(symbols_gpl570), unique(symbols_gpl2700))
	return(length(symbols_alls))
}

# calcula estatisticas 
cutoff.df <- data.frame(cutoff=seq(0, 0.15, length.out=5000))
cutoff.df[, "EarlyDengue"] <- pbapply(cutoff.df, 1, get_intersection_size_adjpv, comp="EarlyDengue vs Control")
cutoff.df[, "LateDengue"] <- pbapply(cutoff.df, 1, get_intersection_size_adjpv, comp="LateDengue vs Control")
cutoff.df[, "Dengue"] <- pbapply(cutoff.df, 1, get_intersection_size_adjpv, comp="Dengue vs Control")

cutoff.melted <- cutoff.df %>% gather("Comparison", "DEGs", c(EarlyDengue, LateDengue, Dengue)) %>% 
	filter(DEGs != 0)

(pl_adjpv <- ggplot(cutoff.melted, aes(x=cutoff, y=DEGs, color=Comparison)) +
	geom_line(size=1) +
	geom_vline(xintercept=0.05, linetype = "longdash", color="gray") +
	geom_vline(xintercept=0.01, linetype = "longdash", color="gray") +
	geom_hline(yintercept=3, linetype = "longdash", color="gray") +
	geom_hline(yintercept=8, linetype = "longdash", color="gray") +
	labs(x="Adjusted P-value", title="Number of Lnoncoding vs Cutoff Used") +
	theme_bw())
for(ext in c("png", "pdf", "svg")){
	ggsave(filename=file.path("figures", "DEG_cutoff", paste0("adjusted_pvalue.", ext)), plot=pl_adjpv)
}


