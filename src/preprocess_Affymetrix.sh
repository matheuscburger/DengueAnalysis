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


echo "Pre processing GSE43777 ... "
src/remove_outliers_Affymetrix.sh GSE43777 GPL570 2>&1 > log/preprocess/remove_outliers_Affymetrix_GSE43777.txt || { echo "Unable to remove outliers for GSE43777"; exit 1; }
echo "Pre processing GSE51808 ... "
src/remove_outliers_Affymetrix.sh GSE51808 GPL13158 2>&1 > log/preprocess/remove_outliers_Affymetrix_GSE51808.txt || { echo "Unable to remove outliers for GSE51808"; exit 1; }
echo "preprocess_Affymetrix done."
