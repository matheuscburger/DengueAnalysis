#!/usr/bin/env python3
# vim:fileencoding=utf8

__author__ = "Matheus Carvalho Bürger"
__email__ = "matheus.cburger@gmail.com"
__license__ = "GPL"

import xlrd
import sys

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit('Usage: %s <detection_pval_file> <expression_file> <sample2gsm_file>' % sys.argv[0])
    det_fname = sys.argv[1]
    signal_fname = sys.argv[2]
    # sample2gsm
    s2gsm_fname = sys.argv[3]
    s2gsm_file = open(s2gsm_fname, "w")
    # pegar arquivo de anotacao
    sample_annot_fname = "config/sample_annotation/with_outliers/GSE13052.tsv"
    with open(sample_annot_fname) as sannot:
        header = sannot.readline().strip().split("\t")
        sample_title_idx = header.index("Sample_title")
        sample_geo_acc_idx = header.index("Sample_geo_accession")
        to_geo = dict()
        s2gsm_file.write("SampleName\tGSM\n")
        # ler arquivo e guardar as colunas Sample_title e Sample_geo_accession
        for line in sannot:
            line = line.strip("\n")
            values = line.split("\t")
            sample_name = values[sample_title_idx].split("-")[-1]
            s2gsm_file.write(sample_name+"\t"+values[sample_geo_acc_idx]+"\n")
            to_geo[sample_name] = values[sample_geo_acc_idx]
    s2gsm_file.close()

    # pegar arquivo bruto GSE13052_illumina_raw.xls
    non_norm_fname = "data/geo_raw/supplemental_data/GSE13052/GSE13052_illumina_raw.xls"
    # abre arquivos de output
    detection_file = open(det_fname, "w")
    signal_file = open(signal_fname, "w")
    tab = xlrd.open_workbook(non_norm_fname)
    ws = tab.sheet_by_index(0)
    num_rows = ws.nrows
    curr_row = 1
    header = [str(v.value) for v in ws.row(0)]
    signal_gsms = []
    pval_gsms = []
    signal_indexes = []
    pval_indexes = []
    for i, h in enumerate(header):
        if h.startswith("AVG_Signal-"):
            sample_name = h.replace("AVG_Signal-", "")
            gsm = to_geo[sample_name]
            signal_gsms.append(gsm)
            signal_indexes.append(i)
        if h.startswith("Detection-"):
            sample_name = h.replace("Detection-", "")
            gsm = to_geo[sample_name]
            pval_gsms.append(gsm)
            pval_indexes.append(i)
    signal_file.write("ProbeName\t"+"\t".join(signal_gsms)+"\n")
    detection_file.write("ProbeName\t"+"\t".join(pval_gsms)+"\n")
    while curr_row < num_rows:
        row = ws.row(curr_row)
        str_row = [str(v.value) for v in row]
        # get and print signal values
        signal_values = [str_row[idx] for idx in signal_indexes]
        signal_values.insert(0, str_row[0])
        signal_file.write("\t".join(signal_values)+"\n")
        # get and print detection p-values
        detection_values = [str_row[idx] for idx in pval_indexes]
        detection_values.insert(0, str_row[0])
        detection_file.write("\t".join(detection_values)+"\n")
        curr_row += 1
    detection_file.close()
    signal_file.close()
