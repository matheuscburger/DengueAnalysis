
library("readr")
library("stringr") 
library("dplyr")
library("corrplot")

dir.create("figures/CEMiTool_joined_corrplot")

fgsea_files <- list.files("results/CEMiTool_joined/fgsea", full.names=T)
nes_files <- grep("NES", fgsea_files, value=T)
pval_files <- grep("p_value", fgsea_files, value=T)


nes_l <- lapply(nes_files, read_tsv)

names(nes_l) <- str_extract(basename(nes_files), "GSE\\d+")

for(n in names(nes_l)){
	nes_l[[n]] <- nes_l[[n]] %>% mutate(Study=n)
}

nes_df <- as.data.frame(do.call(rbind, nes_l))

rownames(nes_df) <- apply(nes_df[, c("pathway", "Study")], 1, paste, collapse="_")

nes_df <- nes_df[order(nes_df[, "pathway"]), setdiff(colnames(nes_df), c("Study", "pathway"))]

nes_mat <- as.matrix(nes_df)



pdf("figures/CEMiTool_joined_corrplot/NES.pdf")
corrColors <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
corrplot(nes_mat, col=corrColors, is.corr=FALSE, addgrid.col="white", insig="blank",
         pch.cex=0.5, pch.col="black", tl.col="black", tl.cex=0.5, cl.cex=0.4, cl.ratio=0.5,
         cl.pos="b", cl.align.text="l", mar=c(0,0,0,0), cl.lim=c(-4,4))
dev.off()
