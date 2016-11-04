#!/usr/bin/env Rscript
library("readr")
library("tidyr")
library("dplyr")
library("biomaRt")

in_args <- commandArgs(trailingOnly=TRUE)
input_fname <- in_args[1]
output_dir <- dirname(input_fname)
gene_col <- in_args[2]

# input_fname <- "/home/mburger/repository/DengueAnalysis/data/processed/CEMiTool_input/GSE28405.tsv"
# output_dir <- "/home/mburger/repository/DengueAnalysis/data/processed/CEMiTool_input/"
# gene_col <- "Symbol"

message("Input file name:", input_fname)
input <- read_tsv(input_fname)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="ensembl.org")
gsymbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "external_gene_name"), filters="ensembl_gene_id",
                  values=input[[gene_col]], mart=ensembl)
by_vector <- "ensembl_gene_id"
names(by_vector) <- gene_col
input <- input %>% left_join(gsymbols, by=by_vector)
res <- input[["hgnc_symbol"]]
res[which(is.na(res) | res == "")] <- input[which(is.na(res) | res == ""), ][["external_gene_name"]]
res[which(is.na(res) | res == "")] <- input[which(is.na(res) | res == ""), ][[gene_col]]
write_tsv(cbind(input[, c(gene_col, "hgnc_symbol", "external_gene_name")], final_name=res), 
		  file.path(output_dir, "gene_conversion", paste0(basename(input_fname))))
input[[gene_col]] <- res
input %>% dplyr::select(-hgnc_symbol, -external_gene_name) %>%
	write_tsv(file.path(output_dir, "annotated", paste0(basename(input_fname))))
