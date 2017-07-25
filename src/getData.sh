#!/bin/bash

for cdir in log data/raw_data data/geo_raw/author data/geo_raw/supplemental_data data/geo_raw/raw_data data/geo_raw/sample_annotation data/geo_raw/probe_annotation; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done


# $? eh o status do ultimo comando
GetGEO() {
	echo "Downloading $1 ..."
	src/microarrayAnalysis/get_geo.R $1 $2 --author-dir=data/geo_raw/author/ --sup-dir=data/geo_raw/supplemental_data/ --sample-annot-dir=data/geo_raw/sample_annotation/ --probe-annot-dir=data/geo_raw/probe_annotation/ --raw-dir=data/geo_raw/raw_data/ 1>&2 2> log/getdata/get_geo_$1.log 
	while [ $? -ne 0 ]; do
		echo "Trying again ..."
		src/microarrayAnalysis/get_geo.R $1 $2 --author-dir=data/geo_raw/author/ --sup-dir=data/geo_raw/supplemental_data/ --sample-annot-dir=data/geo_raw/sample_annotation/ --probe-annot-dir=data/geo_raw/probe_annotation/ --raw-dir=data/geo_raw/raw_data/ 1>&2 2> log/getdata/get_geo_$1.log 
	done	
}


GetGEO GSE43777 GPL570 
GetGEO GSE51808 GPL13158 
GetGEO GSE13052 GPL2700
GetGEO GSE28405 GPL2700

src/getDataGSE28405.py data/raw_data/GSE28405_detection_pval.tsv data/raw_data/GSE28405_signal.tsv data/raw_data/GSE28405_samplename2gsm.tsv || { echo "Unable to extract data for GSE28405"; exit 1; }
src/getDataGSE13052.py data/raw_data/GSE13052_detection_pval.tsv data/raw_data/GSE13052_signal.tsv data/raw_data/GSE13052_samplename2gsm.tsv || { echo "Unable to extract data for GSE13052"; exit 1; }

echo "getData.sh Done."
