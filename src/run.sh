#!/bin/bash

for cdir in log/preprocess log/getdata log/degs log/correlation log/cemitool log/enrichment log/figures; do
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
src/getData.sh &> log/getdata/getData.log

echo $(date +%d-%m-%Y:%H:%M:%S) "Preprocessing data ..."
src/preprocess.sh &> log/preprocess/preprocess.log

echo $(date +%d-%m-%Y:%H:%M:%S) "Running differential expression analysis ..."
src/DEG.sh &> log/degs/DEG.log

echo $(date +%d-%m-%Y:%H:%M:%S) "Running over representation analysis ..."
src/enrichment.sh &> log/enrichment/enrichment.log

echo $(date +%d-%m-%Y:%H:%M:%S) "Running correlation analysis ..."
src/corr.sh &> log/correlation/corr.log

echo $(date +%d-%m-%Y:%H:%M:%S) "Running CEMiTool ..."
src/CEMiTool.sh &> log/cemitool/cemitool.log

echo $(date +%d-%m-%Y:%H:%M:%S) "Getting figures ..."
for script in src/figures/*; do 
	echo $(date +%d-%m-%Y:%H:%M:%S) "Running $script ..."
	base=$(basename $script)
	$script &> log/figures/$base.log
done

echo "Done."
