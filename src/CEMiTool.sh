#!/bin/bash

for cdir in results/CEMiTool data/processed/CEMiTool_input/annotated data/processed/CEMiTool_input/gene_conversion data/processed/CEMiTool_input/filtered tmp/modules results/CEMiTool_joined/enrichment/; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done

# filter using Gustavo's script
parallel -j 10 "src/microarrayAnalysis/genefilter.R --expr-file {} --output {/} --outdir data/processed/CEMiTool_input/filtered --method V --pvalue 0.3 > log/genefilter_{/.}.txt" ::: data/processed/collapsed/*.tsv

# prepare input for CEMiTool
parallel "src/annotateCEMiTool_input.R {} data/processed/CEMiTool_input/ Symbol" ::: data/processed/CEMiTool_input/filtered/GSE*.tsv

#run CEMiTool
parallel -j 1 "src/microarrayAnalysis/CEMiTool.R  {} -o results/CEMiTool/{/.}_CEMiTool -s config/pathways/BTM.gmt -t config/sample_annotation/{/.}.tsv --gene-column Symbol --samples-column Sample_geo_accession -c spearman -n 1000 > log/CEMiTool_{/.}.txt 2> log/CEMiTool_{/.}.txt" :::  data/processed/CEMiTool_input/annotated/*.tsv

echo  "Joining cemitool results ..."
src/microarrayAnalysis/integrateCemitResult.R --output results/CEMiTool_joined.tsv results/CEMiTool/*_GenesInModules.txt

echo "Getting connected components in joined cemitool results ..."
src/microarrayAnalysis/get_modules.py --input results/CEMiTool_joined.tsv --from-col G1 --to-col G2 --filter-col sumOfPairs --filter-val 3 --output results/CEMiTool_joined_modules.txt

echo "Putting each module in one file ..."
parallel 'echo {} | tr "\t" "\n" > tmp/modules/mod{#}.txt' :::: results/CEMiTool_joined_modules.txt

cho "Getting all genesets in enrichr ..."
genesets=$(cat config/libraries_enrichr.txt | sed "s/^/--gs=/" | tr "\n" " ")

echo "Running enrichr ..."
parallel -j 20 "src/microarrayAnalysis/enrichr.py --input {} --output results/CEMiTool_joined/enrichment/{/.}.tsv $genesets" ::: tmp/modules/*.txt

echo "Done."
