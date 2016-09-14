# Informações sobre os estudos utilizando a plataforma GPL2700

Dois estudos que vou usar utilizam a plataforma GPL2700.

Essa plataforma é da empresa Illumina.

Cada estudo libera os dados de uma forma diferente.

## GSE28405
Este estudo possui um arquivo com o nome `GSE28405_non-normalized.txt.gz`.
Neste arquivo a primeira coluna é ID_REF que possui o identificador da probe (e.g. GI_10047089-S).
As próximas duas colunas são "C1" e "Detection Pval".
Suponho que "C1" seja o nome da amostra e "Detection Pval" seja o p-valor de detecção.
Se isso estiver certo existem dois padrões para o nome da amostra: 
- C(?P<n>\d+) - sendo a amostra controle número "n" (1-26)
- P(?P<pat>\d+)T(?P<time>\d+) - sendo a amostra do paciente número "pat" (1-31) e timepoint "time" (1-3)

Para este estudo obter duas tabela com a expressão das amostras e com o "Detection Pval",
mudar o identificador para os identificadores (GSM) do GEO.

## GSE13052
Este estudo possui um arquivo com o nome `GSE13052_illumina_raw.xls`.
Neste arquivo a primeira coluna é TargetID que possui o identificador da probe (e.g. GI_10047089-S).
As próximas colunas são:
- MIN_Signal-1387956034_A
- AVG_Signal-1387956034_A - sinal médio da amostra 1387956034_A
- MAX_Signal-1387956034_A
- NARRAYS-1387956034_A
- ARRAY_STDEV-1387956034_A
- BEAD_STDEV-1387956034_A
- Avg_NBEADS-1387956034_A
- Detection-1387956034_A - p-valor de detecção da amostra 1387956034_A

Para este estudo obter duas tabelas com o AVG_Signal e Detection,
mudar o identificador para os identificadores (GSM) do GEO.
