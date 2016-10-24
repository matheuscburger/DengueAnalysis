#!/bin/bash

for cdir in results/correlation/bygene results/correlation/byprobe; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done

# Calcula correlacao
parallel -j 4 "src/microarrayAnalysis/cor_lnc.R --ovlp config/reannotation/ovlp_lnc.tsv --exp {} --output results/correlation/byprobe/{/} --method spearman 2> log/correlation_{/.}.txt" ::: data/processed/filtered/GSE*.tsv

# Junta correlacao dos 4 estudos
R CMD BATCH src/join_corr.R 
