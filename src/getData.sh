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

GetGEO() {
	echo "Downloading $1 ..."
	src/microarrayAnalysis/get_geo.R $1 $2 --author-dir=data/geo_raw/author/ --sup-dir=data/geo_raw/supplemental_data/ --sample-annot-dir=data/geo_raw/sample_annotation/ --probe-annot-dir=data/geo_raw/probe_annotation/ --raw-dir=data/geo_raw/raw_data/ 1>&2 2> log/get_geo_$1.log 
}


GetGEO GSE43777 GPL570
GetGEO GSE51808 GPL13158
GetGEO GSE13052 GPL2700
GetGEO GSE28405 GPL2700

src/getDataGSE28405.py data/raw_data/GSE28405_detection_pval.tsv data/raw_data/GSE28405_signal.tsv data/raw_data/GSE28405_samplename2gsm.tsv 
src/getDataGSE13052.py data/raw_data/GSE13052_detection_pval.tsv data/raw_data/GSE13052_signal.tsv data/raw_data/GSE13052_samplename2gsm.tsv
