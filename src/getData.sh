#!/bin/bash

for cdir in data/geo_raw/author data/geo_raw/supplemental_data data/geo_raw/raw_data data/geo_raw/sample_annotation data/geo_raw/probe_annotation; do
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
	src/microarrayAnalysis/get_geo.R $1 --author-dir=data/geo_raw/author/ --sup-dir=data/geo_raw/supplemental_data/ --sample-annot-dir=data/geo_raw/sample_annotation/ --probe-annot-dir=data/geo_raw/probe_annotation/ --raw-dir=data/geo_raw/raw_data/ 1>&2 2> log/get_geo_$1.log 
}

export -f GetGEO

parallel -j 10 GetGEO ::: GSE43777 GSE51808 GSE13052 GSE28405

src/getDataGSE28405.py data/raw_data/GSE28405_detection_pval.tsv data/raw_data/GSE28405_signal.tsv data/raw_data/GSE28405_samplename2gsm.tsv 
