#!/bin/bash

for cdir in results/DEG; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done

echo "Obtendo DEGs ... "
parallel -j 10 "src/microarrayAnalysis/do_comparisons.R --comp-file {} --norm-file data/processed/filtered/{/} --stat results/DEG/{/} --sample-annot config/sample_annotation/{/} --annotation-cols Symbol 2> log/do_comparisons_{/.}.txt " ::: config/comparisons/GSE*.tsv

echo "Juntando lista de DEGs ..."
DEG_files=$(for i in results/DEG/GSE*.tsv; do echo "--input $i"; done | tr "\n" " ")
src/microarrayAnalysis/join_stats.R results/joined_DEG.tsv --gene-col Symbol $DEG_files --cores 6
