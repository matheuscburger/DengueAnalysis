
for cdir in figures/enrichment_degs; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done


echo "Enrichment plots for enrichr ..."
#parallel -j 10 "src/microarrayAnalysis/enrichment_plots.R --input {} --output figures/enrichment_degs/enrichr_{/.}.pdf" ::: results/joined_degs/enrichment/enrichr/*.tsv

echo "Enrichment plots for BTMs ..."
parallel -j 10 "src/microarrayAnalysis/enrichment_plots.R --input {} --output figures/enrichment_degs/BTM_{/.}.pdf --term-col ID --pv-col p.adjust --title BTM --rgb-color '#005417' " ::: results/enrichment/do_ora/BTM/*Down*regulated.tsv
parallel -j 10 "src/microarrayAnalysis/enrichment_plots.R --input {} --output figures/enrichment_degs/BTM_{/.}.pdf --term-col ID --pv-col p.adjust --title BTM --rgb-color '#7D0002' " ::: results/enrichment/do_ora/BTM/*Up*regulated.tsv


echo "Enrichment plots for Reactome ..."
parallel -j 10 "src/microarrayAnalysis/enrichment_plots.R --input {} --output figures/enrichment_degs/Reactome_{/.}.pdf --term-col ID --pv-col p.adjust --title Reactome --rgb-color '#005417' " ::: results/enrichment/do_ora/ReactomePathways/*Down*regulated.tsv
parallel -j 10 "src/microarrayAnalysis/enrichment_plots.R --input {} --output figures/enrichment_degs/Reactome_{/.}.pdf --term-col ID --pv-col p.adjust --title Reactome --rgb-color '#7D0002' " ::: results/enrichment/do_ora/ReactomePathways/*Up*regulated.tsv

