#!/bin/bash

for cdir in data/processed/annotated data/processed/collapsed/ data/processed/filtered/; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done

# Annotate
# GPL2700
# pegar somente GPL2700 (Platform == GPL2700)
# pegar somente hits unicos ( Hits == 1 )
# pegar somente anotacao unica (NumAnnot == 1)
# pegar somente coluna de Probe e Gene 
echo "Getting reannotation for GPL2700 ..."
src/get_annotation.py --reannotation config/reannotation/annotation_long.tsv --platform GPL2700 > tmp/reannotation_GPL2700.tsv
echo "Annotating GSE13052 ..."
src/microarrayAnalysis/annotate_probes.R data/processed/normalized/GSE13052.tsv data/processed/annotated/GSE13052.tsv --annotation-file=tmp/reannotation_GPL2700.tsv
echo "Annotating GSE28405 ..."
src/microarrayAnalysis/annotate_probes.R data/processed/normalized/GSE28405.tsv data/processed/annotated/GSE28405.tsv --annotation-file=tmp/reannotation_GPL2700.tsv

# GPL570
echo "Getting reannotation for GPL570 ..."
src/get_annotation.py --reannotation config/reannotation/annotation_long.tsv --platform GPL570 > tmp/reannotation_GPL570.tsv
echo "Annotating GSE43777 ..."
src/microarrayAnalysis/annotate_probes.R data/processed/normalized/GSE43777.tsv data/processed/annotated/GSE43777.tsv --annotation-file=tmp/reannotation_GPL570.tsv
echo "Annotating GSE51808 ..."
src/microarrayAnalysis/annotate_probes.R data/processed/normalized/GSE51808.tsv data/processed/annotated/GSE51808.tsv --annotation-file=tmp/reannotation_GPL570.tsv

## Collapse
#echo "Collapsing ... "
#parallel -j 10 "src/microarrayAnalysis/collapse.R {} data/processed/collapsed/{/} --by-col=Symbol --method=maxmean --annotation-cols=ProbeName" ::: data/processed/annotated/GSE*.tsv

# Filter
echo "Filtering ..."

parallel -j 10 "src/microarrayAnalysis/filter.R {} data/processed/filtered/{/} --method=mean --prop=0.8 --annotation-cols ProbeName --annotation-cols=Symbol" ::: data/processed/annotated/GSE*.tsv

# Collapse
echo "Collapsing ..."

parallel -j 10 "src/microarrayAnalysis/collapse.R {} data/processed/collapsed/{/} --by-col=Symbol --method=maxmean --annotation-cols=ProbeName --annotation-cols=Symbol" ::: data/processed/filtered/GSE*.tsv

echo "Done."
