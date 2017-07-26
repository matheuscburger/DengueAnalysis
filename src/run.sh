#!/bin/bash

for cdir in log/preprocess log/getdata log/degs log/correlation log/cemitool log/enrichment log/figures log/gsea; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done

outliers=''

while [ "$outliers" != Y ] && [ "$outliers" != N ]; do
	echo "Deseja incluir outliers novamente (Y/N) ?"
	read outliers
done

if [ "$outliers" = Y ]; then
	echo "Removing outliers.txt file ..."
	rm config/sample_annotation/outliers.txt
	echo "Copying sample annotation files including outlier samples ..."
	cp config/sample_annotation/with_outliers/* config/sample_annotation/
fi


echo $(date +%d-%m-%Y:%H:%M:%S) " Getting data ..."
src/getData.sh &> log/getdata/getData.log || { echo "Unable to get data"; exit 1; }

echo $(date +%d-%m-%Y:%H:%M:%S) "Preprocessing data ..."
src/preprocess.sh &> log/preprocess/preprocess.log || { echo "Unable to preprocess data"; exit 1; }

echo $(date +%d-%m-%Y:%H:%M:%S) "Running differential expression analysis ..."
src/DEG.sh &> log/degs/DEG.log || { echo "Unable to get DEGs"; exit 1; } 

echo $(date +%d-%m-%Y:%H:%M:%S) "Running over representation analysis ..."
src/enrichment.sh &> log/enrichment/enrichment.log || { echo "Unable to do ORA analysis"; exit 1; }

echo $(date +%d-%m-%Y:%H:%M:%S) "Running gene set enrichment analysis ..."
src/gsea.sh &> log/gsea/gsea.log || { echo "Unable to run GSEA"; exit 1; }

echo $(date +%d-%m-%Y:%H:%M:%S) "Running correlation analysis ..."
src/corr.sh &> log/correlation/corr.log || { echo "Unable to run correlation analysis"; exit 1; }

echo $(date +%d-%m-%Y:%H:%M:%S) "Running CEMiTool ..."
src/CEMiTool.sh &> log/cemitool/cemitool.log || { echo "Unable to run CEMiTool"; exit 1; }

echo $(date +%d-%m-%Y:%H:%M:%S) "Getting figures ..."
for script in src/figures/*; do 
	echo $(date +%d-%m-%Y:%H:%M:%S) "Running $script ..."
	base=$(basename $script)
	$script &> log/figures/$base.log || { echo "Unable to get figures ($script)"; exit 1; }
done

echo $(date +%d-%m-%Y:%H:%M:%S) "Done."
