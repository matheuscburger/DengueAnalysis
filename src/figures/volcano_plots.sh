
for cdir in figures/volcano_plots; do
	echo "Creating $cdir ..."
	if [ -d $cdir ]; then
		echo "$cdir already exists !"
	else
		mkdir -p $cdir
		echo "$cdir created !"
	fi
done

parallel -j 4 "src/microarrayAnalysis/volcano_plots.R --input {} --output figures/volcano_plots/{/.} --y-min=0 --y-max=10 --x-min=-4 --x-max=4" ::: results/DEG/*.tsv
