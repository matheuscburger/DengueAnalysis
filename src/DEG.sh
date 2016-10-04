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
parallel -j 10 "src/microarrayAnalysis/do_comparisons.R --comp-file {} --norm-file data/processed/filtered/{/} --stat results/DEG/{/} --sample-annot config/sample_annotation/{/} --annotation-cols ProbeName --annotation-cols Symbol 2> log/do_comparisons_{/.}.txt " ::: config/comparisons/GSE*.tsv

#echo "Juntando lista de DEGs ..."
#DEG_files=$(for i in results/DEG/GSE*.tsv; do echo "--input $i"; done | tr "\n" " ")
#src/microarrayAnalysis/join_stats.R results/joined_DEG.tsv --gene-col Symbol $DEG_files --cores 6

src/microarrayAnalysis/join_stats.R results/GPL2700_joined_DEG.tsv --probe-col ProbeName --gene-col Symbol --input results/DEG/GSE13052.tsv --input results/DEG/GSE28405.tsv --cores 6 
src/microarrayAnalysis/join_stats.R results/GPL570_joined_DEG.tsv --probe-col ProbeName --gene-col Symbol --input results/DEG/GSE43777.tsv --input results/DEG/GSE51808.tsv --cores 6 
