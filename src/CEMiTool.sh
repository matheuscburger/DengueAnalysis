#!/bin/bash

for cdir in results/CEMiTool tmp/modules tmp/modules_biomart/ results/CEMiTool_joined/enrichment/do_ora/ results/CEMiTool_joined/enrichment/enrichr results/CEMiTool_joined/fgsea figures/enrichment_cemitool_modules/fgsea/ results/CEMiTool_joined/corr_lnc/ results/CEMiTool_joined/corr_stats_lnc/ results/CEMiTool_joined/corr_pvalue_lnc/ results/CEMiTool_joined/corr_padj_lnc/; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done

CEMITOOL=$(Rscript -e "cat(system.file('exec/CEMiTool.R',package='CEMiTool'))")

#run CEMiTool
parallel -j 1 "$CEMITOOL {} --output=results/CEMiTool/{/.}_CEMiTool --sample-annot config/sample_annotation/{/.}.tsv --gene-column Symbol --samples-column Sample_geo_accession --correlation pearson --class-column=Class 2> log/cemitool/CEMiTool_{/.}.txt 1>&2" :::  data/processed/collapsed/GSE*.tsv # || { echo "Unable to run CEMiTool"; exit 1; }

echo  "Joining cemitool results ..."
module_tables=$(find "results/CEMiTool/" -name "module.tsv")
names=$(echo $module_tables | sed "s/\S*\(GSE[0-9]\+\)\S*/--names=\1/g")
input=$(echo $module_tables | sed "s/\(\S\+\)/--input=\1/g")
src/microarrayAnalysis/join_cemitool.R $input $names --output results/CEMiTool_joined.tsv # || { echo "Unable to join CEMiTool"; exit 1; }

echo "Getting connected components in joined cemitool results ..."
#src/microarrayAnalysis/get_modules.py --input results/CEMiTool_joined.tsv --from-col Gene1 --to-col Gene2 --filter-col Sum --filter-val 3 --output results/CEMiTool_joined_modules.txt
src/microarrayAnalysis/get_modules_CEMiTool_joined.R --input results/CEMiTool_joined.tsv --tsv=results/CEMiTool_joined_modules.tsv --graphml=results/CEMiTool_joined_modules.graphml --gmt results/CEMiTool_joined_modules.gmt --chosen-method=walktrap --min-studies=1 --algorithms=fast_greedy --algorithms=label_prop --algorithms=leading_eigen --algorithms=louvain --algorithms=walktrap || { echo "Unable to get modules for CEMiTool joined"; exit 1; }

docker run -ti --rm -e DISPLAY=$DISPLAY -v `pwd`:`pwd` -w `pwd` -v /tmp/.X11-unix:/tmp/.X11-unix tiagopeixoto/graph-tool python3 src/plot_CEMiTool_joined.py results/CEMiTool_joined_modules.graphml results/graphs/clustering_sum_squared_ || echo "Unable to plot graph"

#echo "Generating a GMT file based on connected components ..."
#cat  results/CEMiTool_joined_modules.txt | awk '{ print "Mod"NR"\t\t"$0 }' > results/CEMiTool_joined_modules.gmt

echo "Putting each module in one file ..."
parallel 'echo {} | tr "\t" "\n" > tmp/modules/mod{#}.txt' :::: results/CEMiTool_joined_modules.txt || { echo "Unable to extract modules"; exit 1; }

echo "Converting ensembl ids to hgnc symbol ..."
parallel "src/microarrayAnalysis/ensembl2symbol.R --input {} --output tmp/modules_biomart/{/.}.tsv --is-list" ::: tmp/modules/*.txt || { echo "Unable to get symbols"; exit 1; }

echo "Extracting Gene Symbols ..."
parallel "csvcut -t -c hgnc_symbol {} | sed 1d > tmp/modules_symbols/{/.}.txt" ::: tmp/modules_biomart/*.tsv || { echo "Unable to extract gene symbols"; exit 1; }

echo "Running do_ora for the BTMs ..."
parallel -j 20 "src/microarrayAnalysis/do_ora.R results/CEMiTool_joined/enrichment/do_ora/BTM_{/.}.tsv --genes {} --gmt config/pathways/BTM.gmt" ::: tmp/modules_symbols/*.txt || { echo "Unable to run ORA for CEMiTool joined modules"; exit 1; }

echo "Getting all genesets in enrichr ..."
genesets=$(cat config/libraries_enrichr.txt | sed "s/^/--gs=/" | tr "\n" " ")

echo "Running enrichr ..."
# parallel -j 20 "src/microarrayAnalysis/enrichr.py --input {} --output results/CEMiTool_joined/enrichment/enrichr/{/.}.tsv $genesets" ::: tmp/modules_symbols/*.txt || { echo "Unable to run EnrichR for CemiTool modules"; exit 1; }


echo "Comparisons to GSEA input ..."
parallel -j 10 "src/microarrayAnalysis/comparisons2gseainput.R --dontconvert --input {} --output tmp/fgsea/Log2FC/{/} "  ::: results/DEG/*.tsv || { echo "Unable to get Input for GSEA"; exit 1; }

echo "Running FGSEA ..."
parallel --progress -j 10 "src/microarrayAnalysis/fgsea.R --input {} --output results/CEMiTool_joined/fgsea/{/} --gmt results/CEMiTool_joined_modules.gmt --symbols Symbol" ::: tmp/fgsea/Log2FC/*.tsv || { echo "Unable to run fgsea"; exit 1; }
parallel --progress -j 10 "src/microarrayAnalysis/corrplot_fgsea.R {} figures/enrichment_cemitool_modules/fgsea/{/.}.pdf" ::: results/CEMiTool_joined/fgsea/*.tsv || { echo "Unable to get corrplot"; exit 1; }

echo "Correlation between lncRNAs and genes in modules ..."
parallel -j 10 "src/microarrayAnalysis/lnc_cor_modules.R --modules results/CEMiTool_joined_modules.gmt --expression={} --gencode config/reannotation/gencode_annotation.tsv --biotype config/reannotation/biotypes.tsv --output-corr=results/CEMiTool_joined/corr_lnc/{/} --output-pvalue=results/CEMiTool_joined/corr_pvalue_lnc/{/} --output-padj=results/CEMiTool_joined/corr_padj_lnc/{/} --output-stats=results/CEMiTool_joined/corr_stats_lnc/{/}" ::: data/processed/collapsed/GSE*.tsv || { echo "Unable to calculate correlation between lncRNAs and mRNAs"; exit 1; }

echo "Joining correlation files ..."
src/microarrayAnalysis/join_lnc_cor_modules.R --output results/CEMiTool_joined/corr_lnc_joined.tsv results/CEMiTool_joined/corr_stats_lnc/GSE*.tsv || { echo "Unable to join correlations"; exit 1; }

echo "CEMiTool.sh Done."
