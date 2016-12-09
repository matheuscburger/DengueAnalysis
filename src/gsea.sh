for cdir in tmp/fgsea/annotated tmp/fgsea/colapsed results/fgsea; do
    echo "Creating $cdir ..."
    if [ -d $cdir ]; then
        echo "$cdir already exists !"
    else
        mkdir -p $cdir
        echo "$cdir created !"
    fi
done

parallel -j 10 "src/microarrayAnalysis/ensembl2symbol.R --input {} --output tmp/fgsea/annotated/{/} --ensembl-col Symbol"  ::: data/processed/filtered/*.tsv

parallel -j 10 "src/microarrayAnalysis/collapse.R {} tmp/fgsea/colapsed/{/} --by-col external_gene_name --method maxmean --annotation-cols  hgnc_symbol --annotation-cols ensembl_gene_id --annotation-cols ProbeName" ::: tmp/fgsea/annotated/*.tsv

for database in config/pathways/*.gmt; do
    parallel --progress -j 10 "mkdir -p results/fgsea/{/.}/$database; src/microarrayAnalysis/fgsea.R --expression {} --es results/fgsea/{/.}/$database/es.tsv --nes results/fgsea/{/.}/$database/nes.tsv --pval results/fgsea/{/.}/$database/pval.tsv --gmt config/pathways/ReactomePathways.gmt --sample-annotation config/sample_annotation/{/} --symbols external_gene_name --sample-name-col Sample_geo_accession --class-col Class" ::: tmp/fgsea/colapsed/*.tsv
done
