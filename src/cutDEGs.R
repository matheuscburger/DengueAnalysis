
library("readr")
library("dplyr")
library("tidyr")
library("stringr")
library("biomaRt")

# le GPL570
gpl570 <- read_tsv("results/GPL570_joined_DEG.tsv", guess_max=100000)
gpl2700 <- read_tsv("results/GPL2700_joined_DEG.tsv", guess_max=10000000)

# conta o numero de estudos por teste
countStudies <- function(x){
	# pega Log2FC
	log2fccols <- grep("Log2FC", colnames(x), value=T)
	countl2fc <- apply(apply(x[, log2fccols], 2, is.na), 1, function(x) sum(!x))
	# ve quantos NAs por linha
	x[, "Log2FC_Count"] <- countl2fc
	# pega o maximo para uma determinada comparacao
	nstudies <- sapply(split(countl2fc, x[, "title"]), max)
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


getDEGs <- function(x, padj_cut=1, p_cut=1, lfc_cut=0){
	# Pega colunas para filtrar
	padj_cols <- grep("LIMMA\\.adjp\\.", colnames(x), value=TRUE)
	names(padj_cols) <- str_match(padj_cols, "\\.(\\w*)$")[,2]
	p_cols <- grep("LIMMA\\.rawp", colnames(x), value=TRUE)
	names(p_cols) <- str_match(p_cols, "\\.(\\w*)$")[,2]
	lfc_cols <- grep("Log2FC\\.", colnames(x), value=TRUE)
	names(lfc_cols) <- str_match(lfc_cols, "\\.(\\w*)$")[,2]

	studies <- unique(c(names(lfc_cols), names(p_cols), names(padj_cols)))

	# filters
	padj_filter <- x[, padj_cols, drop=FALSE] <= padj_cut
	colnames(padj_filter) <- names(padj_cols)
	p_filter <- x[, p_cols, drop=FALSE] <= p_cut
	colnames(p_filter) <- names(p_cols)
	lfc_up_filter <- x[, lfc_cols, drop=FALSE] >= lfc_cut
	lfc_down_filter <- x[, lfc_cols, drop=FALSE] <= -lfc_cut
	colnames(lfc_up_filter) <- colnames(lfc_down_filter) <- names(lfc_cols)

	filter_up <- padj_filter[, studies] & p_filter[, studies] & lfc_up_filter[, studies]
	filter_down <- padj_filter[, studies] & p_filter[, studies] & lfc_down_filter[, studies]
	#x[which(filter_boolean), ]
	res_mat <- filter_up
	res_mat[which(!filter_up, arr.ind=T)] <- "Equal"
	res_mat[which(filter_up, arr.ind=T)] <- "Up-regulated"
	res_mat[which(filter_down, arr.ind=T)] <- "Down-regulated"

	res <- x %>% select(title, ProbeName, Symbol) %>% bind_cols(data.frame(res_mat))
	return(res) 
}

collapse_DEGs <- function(x){
	x <- unique(x)
	if(any(x %in% c("Up-regulated", "Down-regulated"))){
		x <- x[x != "Equal"]
	}
	x <- x[!is.na(x)]
	res <- paste0(x, collapse="/")
	return(res)
}

degsGPL2700 <- getDEGs(gpl2700, padj_cut=0.05) %>% group_by(title, Symbol) %>% 
	summarise(GSE13052=collapse_DEGs(GSE13052),
	          GSE28405=collapse_DEGs(GSE28405))

degsGPL570 <- getDEGs(gpl570, padj_cut=0.05) %>% group_by(title, Symbol) %>% 
	summarise(GSE43777=collapse_DEGs(GSE43777),
	          GSE51808=collapse_DEGs(GSE51808))

degs <- full_join(degsGPL570, degsGPL2700)

# conta o numero de estudos por teste
gsecols <- grep("GSE\\d+", colnames(degs), value=T)
countgse <- apply(degs[, gsecols], 1, function(x) sum(!(is.na(x) | x == "")) )
# pega o maximo para uma determinada comparacao
nstudies <- sapply(split(countgse, degs[, "title"]), max)
degs[, "NStudies"] <- nstudies[degs[["title"]]]
degs[, "HasNAs"] <- degs[, "NStudies"] != countgse
degs[, "HasDiscordantDEGs"] <- apply(degs[, gsecols], 1, function(x){any(grepl("/", x))})
degs[, "Direction"] <- apply(degs[, gsecols], 1, function(x){
		  res <- "NotSignificant"
		  down <- sum(x == "Down-regulated", na.rm=T)
		  up <- sum(x == "Up-regulated", na.rm=T)
		  if(down > 0 & up > 0){
			  res <- "Discordant"
		  } else if (down > 0){
			  res <- "Down-regulated"
		  } else if (up > 0){
			  res <- "Up-regulated"
		  }
		  return(res)
			  })

degs[, "CountDEGs"] <- apply(degs[, gsecols], 1, function(x){
		  down <- sum(x == "Down-regulated", na.rm=T)
		  up <- sum(x == "Up-regulated", na.rm=T)
		  return(max(down, up))
			  })

degs <- degs %>% mutate(Ratio=CountDEGs/NStudies)

degs %>% filter(Ratio >= 0.6, !HasDiscordantDEGs) %>%
	count(title, Direction, NStudies) %>%
	write_tsv("results/DEG_stats.tsv")


# Anotar tabela de degs
# le anotacao
reannotation <- read_tsv("config/reannotation/annotation_long.tsv") %>% 
	select(Database, Gene, Type, Group) %>% unique
degs <- degs %>% left_join(reannotation, by=c("Symbol" = "Gene"))

# Stats by group
degs %>% filter(Ratio == 1, !HasDiscordantDEGs) %>%
	count(title, Direction, NStudies, Group) %>%
	write_tsv("results/DEG_stats_by_group.tsv")

ensgs <- degs %>% ungroup() %>% filter(Database == "gencode") %>% dplyr::select(Symbol) %>% unique %>% .[["Symbol"]]

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="ensembl.org")
gsymbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "external_gene_name"), filters="ensembl_gene_id",
                  values=ensgs, mart=ensembl)

degs <- degs %>% left_join(gsymbols, by=c("Symbol"="ensembl_gene_id"))

write_tsv(degs, "results/joined_DEG.tsv")
