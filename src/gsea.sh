for cdir in tmp/fgsea/Log2FC tmp/fgsea/annotated tmp/fgsea/colapsed tmp/fgsea/scaled tmp/fgsea/meanzscore results/fgsea; do
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
	for cdir in results/fgsea/Log2FC/$db results/fgsea/MeanZScore/$db figures/GSEA/Log2FC/$db figures/GSEA/MeanZScore/$db; do
		echo "Creating $cdir ..."
		if [ -d $cdir ]; then
			echo "$cdir already exists !"
		else
			mkdir -p $cdir
			echo "$cdir created !"
		fi
	done
done

echo "Comparisons to GSEA input ..."
parallel -j 10 "src/microarrayAnalysis/comparisons2gseainput.R --input {} --output tmp/fgsea/Log2FC/{/} "  ::: results/DEG/*.tsv

echo "Ensembl to Gene Symbols  ..."
parallel -j 10 "src/microarrayAnalysis/ensembl2symbol.R --input {} --output tmp/fgsea/annotated/{/} --ensembl-col Symbol"  ::: data/processed/filtered/*.tsv

echo "Collapsing by external_gene_name ..."
parallel -j 10 "src/microarrayAnalysis/collapse.R {} tmp/fgsea/colapsed/{/} --by-col external_gene_name --method maxmean --annotation-cols  hgnc_symbol --annotation-cols ensembl_gene_id --annotation-cols ProbeName" ::: tmp/fgsea/annotated/*.tsv

echo "Scaling data (row z-score) ..."
parallel -j 10 "src/microarrayAnalysis/scale.R --expr-file {} --output tmp/fgsea/scaled/{/} --annotation-cols external_gene_name" ::: tmp/fgsea/colapsed/*.tsv

echo "Calculating group averages ..."
parallel -j 10 "src/microarrayAnalysis/group_mean.R --expr-file {} --group-file config/sample_annotation/{/} --output tmp/fgsea/meanzscore/{/} --group-col Stage --samplename-col Sample_geo_accession" ::: tmp/fgsea/scaled/*.tsv

echo "Running FGSEA ..."
for database in config/pathways/*.gmt; do
    db=$(basename $database .gmt)
    parallel --progress -j 10 "src/microarrayAnalysis/fgsea.R --input {} --output results/fgsea/Log2FC/$db/{/} --gmt $database --symbols hgnc_symbol" ::: tmp/fgsea/Log2FC/*.tsv
    parallel --progress -j 10 "src/microarrayAnalysis/fgsea.R --input {} --output results/fgsea/MeanZScore/$db/{/} --gmt $database --symbols external_gene_name" ::: tmp/fgsea/meanzscore/*.tsv

    parallel --progress -j 10 "src/microarrayAnalysis/corrplot_fgsea.R {} figures/GSEA/Log2FC/$db/{/.}.pdf" ::: results/fgsea/Log2FC/$db/*.tsv
    parallel --progress -j 10 "src/microarrayAnalysis/corrplot_fgsea.R {} figures/GSEA/MeanZScore/$db/{/.}.pdf" ::: results/fgsea/MeanZScore/$db/*.tsv
done
