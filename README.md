# [Informação sobre os estudos utilizados](StudyInfo.md)
# Estrutura de diretórios

- / - Raiz do projeto
- /src - Deve (idealmente) conter todas os scripts do projeto
- /src/figures - scripts para gerar as figuras da apresentação do power-point "To-do list"
- /src/microarrayAnalysis -  scripts do pipeline original e outros para fazer a análise de microarray
- /config - Deve conter todas as informações necessárias para gerar todas as outras coisas
- /config/sample_annotation - informação sobre as amostras sem outliers
- /config/sample_annotation/with_outliers - informação sobre as amostras contendo outliers
- /figures - figuras da apresentação do power-point "To-do list"
- /data - contém os dados baixados e pré-processados (não resultados)
- /data/geo_raw - dados baixados do GEO
- /data/geo_raw/author - matriz de expressão gênica pré-processado (ou não) pelo autor
- /data/geo_raw/probe_annotation - anotação das probes proveniente do GEO
- /data/geo_raw/raw_data - arquivos brutos de cada amostra como CEL files
- /data/geo_raw/sample_annotation - arquivos de anotação
- /data/geo_raw/supplemental_data - arquivos suplementares
- /data/raw_data - raw data que veio do GEO mas foi necessário algum parsing, descompactação, etc de arquivos
- /data/processed - dados processados
- /data/processed/normalized - dados normalizados
- /quality_control - para guardar relatórios de controle de qualidade
- /quality_control/before_norm - controle de qualitade antes de normalizar
- /quality_control/after_norm - controle de qualitade depois de normalizar

