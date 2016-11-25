#!/bin/bash

for cdir in results/joined_degs/enrichment/do_ora results/joined_degs/enrichment/enrichr results/joined_degs/genes; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done


csvgrep -t -i -c HasDiscordantDEGs -m "TRUE" results/joined_DEG.tsv |  # get only concordant DEGs
	csvgrep -i -c Direction -m "Discordant" |  # get only concordant DEGs
	awk -F"," '$12 > 0.6' |  # keep only rows with ratio above 0.6 / Degs in at least 3 comparisons
	csvcut -c title,Direction,external_gene_name |  # get only three columns
	awk -F, '{print $1"_"$2"\t"$3}' |  # join first two columns
	sed "s/ /_/g" |  # remove blank spaces
	sed 1d |  # remove first line
	awk '{print $2 > "tmp/joined_degs/"$1".txt"}'  # put each comparison in one file

parallel "sort {} | uniq > results/joined_degs/genes/{/}" ::: tmp/joined_degs/*.txt

echo "Running do_ora for the BTMs ..."
parallel -j 20 "src/microarrayAnalysis/do_ora.R results/joined_degs/enrichment/do_ora/BTM_{/.}.tsv --genes {} --gmt config/pathways/BTM.gmt" ::: results/joined_degs/genes/*.txt

echo "Getting all genesets in enrichr ..."
genesets=$(cat config/libraries_enrichr.txt | sed "s/^/--gs=/" | tr "\n" " ")

echo "Running enrichr ..."
parallel -j 20 "src/microarrayAnalysis/enrichr.py --input {} --output results/joined_degs/enrichment/enrichr/{/.}.tsv $genesets" ::: results/joined_degs/genes/*.txt

echo "Done."
