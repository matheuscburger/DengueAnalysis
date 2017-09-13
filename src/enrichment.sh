#!/bin/bash

for cdir in results/enrichment/do_ora results/enrichment/enrichr results/joined_degs/genes tmp/joined_degs; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done
for database in config/pathways/*.gmt; do
	db=$(basename $database .gmt)
	cdir=results/enrichment/do_ora/$db
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
	awk '{print $2 > "tmp/joined_degs/"$1".txt"}' || { echo "Unable to extract DEGs"; exit 1; } # put each comparison in one file

parallel "sort {} | uniq > results/joined_degs/genes/{/}" ::: tmp/joined_degs/*.txt || { echo "Unable to sort and write DEG list"; exit 1; }

for database in config/pathways/*.gmt; do
    db=$(basename $database .gmt)
	echo "Running do_ora for the BTMs ..."
	parallel -j 20 "src/microarrayAnalysis/do_ora.R results/enrichment/do_ora/$db/{/.}.tsv --genes {} --gmt ${database}" ::: results/joined_degs/genes/*.txt || { echo "Unable to run ORA for $db"; exit 1; }
done


#echo "Getting all genesets in enrichr ..."
#genesets=$(cat config/libraries_enrichr.txt | sed "s/^/--gs=/" | tr "\n" " ")

#echo "Running enrichr ..."
#parallel -j 20 "src/microarrayAnalysis/enrichr.py --input {} --output results/enrichment/enrichr/{/.}.tsv $genesets" ::: results/joined_degs/genes/*.txt ||  { echo "Unable to run Enrichr"; exit 1; }

echo "Done."
