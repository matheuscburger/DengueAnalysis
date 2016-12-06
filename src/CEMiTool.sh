#!/bin/bash

for cdir in results/CEMiTool data/processed/CEMiTool_input/annotated data/processed/CEMiTool_input/gene_conversion data/processed/CEMiTool_input/filtered data/processed/CEMiTool_FGSEA_input/gene_conversion data/processed/CEMiTool_FGSEA_input/annotated tmp/modules results/CEMiTool_joined/enrichment/enrichr results/CEMiTool_joined/enrichment/do_ora results/CEMiTool_joined/fgsea/ tmp/with_symbols; do
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

echo "Generating a GMT file based on connected components ..."
cat  results/CEMiTool_joined_modules.txt | awk '{ print "Mod"NR"\t\t"$0 }' > results/CEMiTool_joined_modules.gmt

echo "Putting each module in one file ..."
parallel 'echo {} | tr "\t" "\n" > tmp/modules/mod{#}.txt' :::: results/CEMiTool_joined_modules.txt


echo "Running do_ora for the BTMs ..."
parallel -j 20 "src/microarrayAnalysis/do_ora.R results/CEMiTool_joined/enrichment/do_ora/BTM_{/.}.tsv --genes {} --gmt config/pathways/BTM.gmt" ::: tmp/modules/*.txt

echo "Getting all genesets in enrichr ..."
genesets=$(cat config/libraries_enrichr.txt | sed "s/^/--gs=/" | tr "\n" " ")

echo "Running enrichr ..."
parallel -j 20 "src/microarrayAnalysis/enrichr.py --input {} --output results/CEMiTool_joined/enrichment/enrichr/{/.}.tsv $genesets" ::: tmp/modules/*.txt

echo "Preparing input for FGSEA ..."
parallel "src/annotateCEMiTool_input.R {} data/processed/CEMiTool_FGSEA_input/ Symbol" ::: data/processed/collapsed/GSE*.tsv

echo "Running FGSEA ..."
parallel -j 10 "src/microarrayAnalysis/fgsea.R --expression {} --es results/CEMiTool_joined/fgsea/ES_{/} --nes results/CEMiTool_joined/fgsea/NES_{/} --pval results/CEMiTool_joined/fgsea/p_value_{/} --gmt results/CEMiTool_joined_modules.gmt --sample-annotation config/sample_annotation/{/} --annotation-cols Symbol --symbols Symbol --sample-name-col Sample_geo_accession --class-col Class" ::: data/processed/CEMiTool_FGSEA_input/annotated/GSE*.tsv

echo "CEMiTool.sh Done."
