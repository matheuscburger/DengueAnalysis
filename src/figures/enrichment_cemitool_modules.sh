
for cdir in figures/enrichment_cemitool_modules; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done


echo "Enrichment plots for enrichr ..."
parallel -j 10 "src/microarrayAnalysis/enrichment_plots.R --input {} --output figures/enrichment_cemitool_modules/enrichr_{/.}.pdf" ::: results/CEMiTool_joined/enrichment/enrichr/*.tsv

echo "Enrichment plots for BTMs ..."
parallel -j 10 "src/microarrayAnalysis/enrichment_plots.R --input {} --output figures/enrichment_cemitool_modules/BTM_{/.}.pdf --term-col ID --pv-col p.adjust --title BTM" ::: results/CEMiTool_joined/enrichment/do_ora/*.tsv
