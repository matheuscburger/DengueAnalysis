# Estrutura de diretórios

- / - Raiz do projeto
- /src - Deve (idealmente) conter todas os scripts do projeto
- /src/figures - scripts para gerar as figuras da apresentação do power-point "To-do list"
- /src/microarrayAnalysis -  scripts do pipeline original e outros para fazer a análise de microarray
- /config - Deve conter todas as informações necessárias para gerar todas as outras coisas
- /config/sample_annotation - informação sobre as amostras

# Coisas feitas
- Criei o diretório `config/sample_annotation` para guardar a informação das amostras
``` 
mkdir -p config/sample_annotation
```

- Criei arquivo `config/studies.tsv` com o ID dos estudos

- Como baixar sample_annotation:
```
parallel "src/microarrayAnalysis/get_sample_annot.py --out config/sample_annotation/{}.tsv {}" ::: $(cat config/studies.tsv | cut -f 1 | sed 1d)
```

- Adicionei manualmente as colunas Class, ExtendedClass e Time (caso estudo for timecourse)
