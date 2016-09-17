if [ ! $# -eq 1 ]; then
	echo "Please give me the study id !"
	exit 1
fi

study=$1

echo "Study $study ..."

# Cria diretorio tmp
if [ ! -d tmp ]; then
	echo "Creating tmp ..."
	mkdir tmp
fi

# Cria arquivo de outliers
if [ ! -f config/sample_annotation/outliers.txt ]; then
	echo "Criando arquivo de outliers ..."
	touch config/sample_annotation/outliers.txt
fi

# dados da Illumina estao no diretorio data/raw_data


cp data/raw_data/${study}_signal.tsv tmp/${study}.tsv

num_outliers=1
count=0
# remove outliers antes da normalizacao
echo "Avaliando outliers antes da normalização ..."
while [ $num_outliers -gt 0 ]; do
	echo "Round number: " $count
	# arrayQualityMetrics antes de normalizar
	src/microarrayAnalysis/aqm.R --input-exp tmp/${study}.tsv config/sample_annotation/${study}.tsv quality_control/before_norm/${study}_$count --int-cols ExtendedClass --int-cols Class 2> log/aqm_before_norm_${study}_$count.txt
	# obtem outliers em pelo menos 2 metodos
	src/microarrayAnalysis/get_outliers_from_json.py quality_control/before_norm/${study}_$count/outliers.json --count=2 > quality_control/before_norm/${study}_$count/outliers.txt 
	cat quality_control/before_norm/${study}_$count/outliers.txt config/sample_annotation/outliers.txt > tmp/outliers.txt
	cat tmp/outliers.txt  | sort | uniq | sed '/^\s*$/d' >> config/sample_annotation/outliers.txt

	# normaliza
	src/microarrayAnalysis/do_quantile.R tmp/${study}.tsv data/processed/normalized/${study}.tsv --annotation-cols=ProbeName 2> log/quantile_${study}_$count.txt

	# arrayQualityMetrics depois de normalizar
	src/microarrayAnalysis/aqm.R --input-exp data/processed/normalized/${study}.tsv config/sample_annotation/${study}.tsv quality_control/after_norm/${study}_$count --int-cols ExtendedClass --int-cols Class 2> log/aqm_after_norm_${study}_$count.txt

	# remove outliers da anotacao das amostras
	parallel "grep -vf config/sample_annotation/outliers.txt {} > config/sample_annotation/{/}" ::: config/sample_annotation/with_outliers/GSE*.tsv
	# remove outliers da tabela de expressao
	cols_to_remove=$(cat config/sample_annotation/outliers.txt | sed "s/^/--columns=/" | tr "\n" " ")
	src/microarrayAnalysis/remove_columns.py data/raw_data/${study}_signal.tsv  $cols_to_remove > tmp/${study}.tsv

	# conta outliers
	num_outliers=$(wc -l quality_control/before_norm/${study}_$count/outliers.txt | sed 's/\s.\+//')
	echo $num_outliers " outliers ..."

	# atualiza contador
	count=$((count+1))
done

# remove outliers depois da normalização
num_outliers=1
echo "Avaliando outliers depois da normalização ..."
while [ $num_outliers -gt 0 ]; do
	echo "Round number: " $count
	# arrayQualityMetrics antes de normalizar
	src/microarrayAnalysis/aqm.R --input-exp tmp/${study}.tsv config/sample_annotation/${study}.tsv quality_control/before_norm/${study}_$count --int-cols ExtendedClass --int-cols Class 2> log/aqm_before_norm_${study}_$count.txt

	# normaliza
	src/microarrayAnalysis/do_quantile.R tmp/${study}.tsv data/processed/normalized/${study}.tsv --annotation-cols=ProbeName 2> log/quantile_${study}_$count.txt

	# arrayQualityMetrics depois de normalizar
	src/microarrayAnalysis/aqm.R --input-exp data/processed/normalized/${study}.tsv config/sample_annotation/${study}.tsv quality_control/after_norm/${study}_$count --int-cols ExtendedClass --int-cols Class 2> log/aqm_before_norm_${study}_$count.txt

	# obtem outliers em pelo menos 2 metodos
	src/microarrayAnalysis/get_outliers_from_json.py quality_control/after_norm/${study}_$count/outliers.json --count=2 > quality_control/after_norm/${study}_$count/outliers.txt 
	cat quality_control/after_norm/${study}_$count/outliers.txt config/sample_annotation/outliers.txt > tmp/outliers.txt
	cat tmp/outliers.txt | sort | uniq | sed '/^\s*$/d' >> config/sample_annotation/outliers.txt

	# remove outliers da anotacao das amostras
	parallel "grep -vf config/sample_annotation/outliers.txt {} > config/sample_annotation/{/}" ::: config/sample_annotation/with_outliers/GSE*.tsv
	# remove outliers da tabela de expressao
	cols_to_remove=$(cat config/sample_annotation/outliers.txt | sed "s/^/--columns=/" | tr "\n" " ")
	src/microarrayAnalysis/remove_columns.py data/raw_data/${study}_signal.tsv  $cols_to_remove > tmp/${study}.tsv

	# conta outliers
	num_outliers=$(wc -l quality_control/after_norm/${study}_$count/outliers.txt | sed 's/\s.\+//')
	echo $num_outliers " outliers ..."

	# atualiza contador
	count=$((count+1))
done
