#!/bin/bash

for cdir in results/CEMiTool data/processed/CEMiTool_input/annotated data/processed/CEMiTool_input/gene_conversion; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done

parallel -j 10 "src/microarrayAnalysis/genefilter.R --expr-file {} --output {/} --outdir data/processed/CEMiTool_input --method V --pvalue 0.3 > log/genefilter_{/.}.txt" ::: data/processed/collapsed/*.tsv
parallel "src/annotateCEMiTool_input.R {} Symbol" ::: data/processed/CEMiTool_input/GSE*.tsv
parallel -j 10 "src/microarrayAnalysis/CEMiTool.R  {} -o results/CEMiTool/{/.}_CEMiTool -s config/pathways/BTM.gmt -t config/sample_annotation/{/.}.tsv --gene-column Symbol --samples-column Sample_geo_accession -c spearman -n 1000 > log/CEMiTool_{/.}.txt 2> log/CEMiTool_{/.}.txt" :::  data/processed/CEMiTool_input/annotated/*.tsv
