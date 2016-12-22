
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
parallel -j 10 "src/microarrayAnalysis/enrichment_plots.R --input {} --output figures/enrichment_degs/BTM_{/.}.pdf --term-col ID --pv-col p.adjust --title BTM --rgb-color '#005417' " ::: results/joined_degs/enrichment/do_ora/*Down*regulated.tsv
parallel -j 10 "src/microarrayAnalysis/enrichment_plots.R --input {} --output figures/enrichment_degs/BTM_{/.}.pdf --term-col ID --pv-col p.adjust --title BTM --rgb-color '#7D0002' " ::: results/joined_degs/enrichment/do_ora/*Up*regulated.tsv
