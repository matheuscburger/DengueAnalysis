#!/usr/bin/env Rscript

library("dplyr")
library("tidyr")
library("readr")
library("scales")
library("ggplot2")
library("gplots")
library("stringr")
library("extrafont")
font_import(pattern="DejaVuSans*", prompt=FALSE)
loadfonts(device="pdf")

dir.create(file.path("figures", "correlation"))
dir.create(file.path("figures", "correlation", "scatter_plots"))
dir.create(file.path("figures", "correlation", "scatter_plots_ranks"))

# load correlations
corr <- read_tsv("results/correlation/all.tsv")

# filtered table
(corr.filtered <- corr %>% 
 	filter(abs(min.cor) > 0.25 & max.adj.p < 0.1))


# auxiliar function to discover the position of lncRNA in relation to mRNA
# preference: overlap > upstream > downstream
aux_fun <- function(x){
	if("overlap" %in% x) {
		return("overlap")
	} else if("upstream" %in% x) {
		return("upstream")
	} else {
		return("downstream")
	}
}
corr.uniq <- corr.filtered %>% group_by(mrn_id, lnc_id) %>% 
	mutate(selected=aux_fun(int_type)) %>% 
	filter(selected == int_type) %>%
	ungroup()

aux <-  corr.uniq %>%
	dplyr::select(starts_with("correlation")) %>% as.matrix

rownames(aux) <- apply(corr.uniq[, c("mrn_id", "lnc_id")], 1, paste0, collapse="_")
#		select(selected, mrn_id, lnc_id, int_type)
pdf("figures/correlation/corr_heatmap.pdf")
heatmap.2(aux, trace="none", col=bluered, margins=c(9, 20))
dev.off()
write.table(aux, file.path("figures", "correlation", "corr_heatmap.aux.tsv"),
            sep="\t", quote=FALSE, col.names=NA) 

# Read expression files
studies <- unique(str_extract(colnames(corr.uniq), "GSE\\d+"))
studies <- studies[!is.na(studies)]
exp_list <- lapply(file.path("data/processed/filtered", paste0(studies, ".tsv")), read_tsv)
names(exp_list) <- studies

# read Sample annotations 
samp_annot_list <- lapply(file.path("config/sample_annotation/", paste0(studies, ".tsv")), read_tsv)
names(samp_annot_list) <- studies

dir.create(file.path("figures", "correlation", "scatter_plots_data"))

# make plots
for(i in 1:nrow(corr.uniq)){
	lnc.name <- corr.uniq[["lnc_id"]][i]
	mrna.name <- corr.uniq[["mrn_id"]][i] 
	name <- paste0(lnc.name, ".vs.", mrna.name)
	rank.plots <- list()
	scatter.plots <- list()
	for(s in studies){
		classes <- samp_annot_list[[s]] %>% dplyr::select(Sample_geo_accession, Class)
		lnc.probename <- corr.uniq[[paste0("ProbeName_lnc.", s)]][i]
		mrna.probename <- corr.uniq[[paste0("ProbeName_mrn.", s)]][i]
		corr_value <- corr.uniq[[paste0("correlation.", s)]][i]
		lnc <- exp_list[[s]] %>% 
			filter(ProbeName == lnc.probename) %>% 
			dplyr::select(-Symbol, -ProbeName) %>% 
			gather("Sample", "lncRNA") 
		mrna <- exp_list[[s]] %>% 
			filter(ProbeName == mrna.probename) %>% 
			dplyr::select(-Symbol, -ProbeName) %>% 
			gather("Sample", "mRNA")
		exp_curr <- inner_join(lnc, mrna, by=c("Sample"="Sample")) %>%
			left_join(classes, by=c("Sample"="Sample_geo_accession"))
		pl <- ggplot(exp_curr, aes(x=rank(lncRNA), y=rank(mRNA))) + 
			geom_point(aes(color=Class), size=3) +
			geom_smooth(method="lm", col="darkblue", se=FALSE) +
			scale_color_manual(values=c("Dengue"="#9F2110", "Control"="#164D9F")) +
			annotate("text", x=Inf, y=-Inf, label=paste("Spearman Correlation = ", round(corr_value, digits=3)), 
					  hjust=1.1, vjust=-1, size=6, family="DejaVu Sans") +
			labs(title=s, x=paste0("rank(", lnc.probename, ")"), y=paste0("rank(", mrna.probename, ")")) +
			theme_bw() +
			theme(text=element_text(size=15, family="DejaVu Sans"))
		rank.plots[[s]] <- pl
		pl <- ggplot(exp_curr, aes(x=lncRNA, y=mRNA)) + 
			geom_point(aes(color=Class), size=3) +
			geom_smooth(method="lm", col="darkblue", se=FALSE) +
			scale_color_manual(values=c("Dengue"="#9F2110", "Control"="#164D9F")) +
			annotate("text", x=Inf, y=-Inf, label=paste("Spearman Correlation = ", round(corr_value, digits=3)), 
					  hjust=1.1, vjust=-1, size=6, family="DejaVu Sans") +
			labs(title=s, x=paste0(lnc.probename), y=paste0(mrna.probename)) +
			theme_bw() +
			theme(text=element_text(size=15, family="DejaVu Sans"))
		scatter.plots[[s]] <- pl
        write_tsv(exp_curr, file.path("figures", "correlation", 
                                      "scatter_plots_data", 
                                      paste0(s, "_", mrna.name, "_", lnc.name, ".tsv")))
	}
	pdf(file.path("figures", "correlation", "scatter_plots_ranks", paste0(name, ".pdf")))
	.null <- sapply(rank.plots, print)
	dev.off()
	pdf(file.path("figures", "correlation", "scatter_plots", paste0(name, ".pdf")))
	.null <- sapply(scatter.plots, print)
	dev.off()

}
