#!/bin/bash

for cdir in data/processed/normalized quality_control/after_norm quality_control/before_norm; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done


echo "Pre processing GSE13052 ... "
src/remove_outliers_Illumina.sh GSE13052 2>&1 > log/preprocess/remove_outliers_Illumina_GSE13052.txt || echo "Unable to remove outliers for GSE23052"
echo "Pre processing GSE28045 ... "
src/remove_outliers_Illumina.sh GSE28405 2>&1 > log/preprocess/remove_outliers_Illumina_GSE28405.txt || echo "Unable to remove outliers for GSE28405"
echo "preprocess_Illumina done."
