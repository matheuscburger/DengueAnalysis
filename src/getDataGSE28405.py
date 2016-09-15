#!/usr/bin/env python3
# vim:fileencoding=utf8

__author__ = "Matheus Carvalho Bürger"
__email__ = "matheus.cburger@gmail.com"
__license__ = "GPL"

import re
import gzip
import sys

if __name__ == "__main__":
    # parametros
    if len(sys.argv) < 4:
            sys.exit('Usage: %s <detection_pval_file> <expression_file> <sample2gsm_file>' % sys.argv[0])
    # det_fname = "detection_pvalue.tsv"
    det_fname = sys.argv[1]
    # signal_fname = "expression.tsv"
    signal_fname = sys.argv[2]
    # sample2gsm
    s2gsm_fname = sys.argv[3]
    s2gsm_file = open(s2gsm_fname, "w")
    # pegar arquivo de anotacao
    sample_annot_fname = "config/sample_annotation/with_outliers/GSE28405.tsv"
    control_re = re.compile("Control_sample\s+(?P<n>\d+)")
    patient_re = re.compile("Patient_Timepoint(?P<time>\d+)\s+(?P<pat>\d+)")
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
            # Pegar Sample_name que sera chave e Sample_geo_accession sera valor
            # Se Control_sample\s+(?P<n>\d+)
            mtch = control_re.search(values[sample_title_idx])
            if mtch:
                # Sample_name sera C\n
                sample_name = "C"+mtch.group("n")
            else:
                # Senão - .*Patient_Timepoint(?P<time>\d+)\s+(?P<pat>\d+)
                mtch = patient_re.search(values[sample_title_idx])
                if mtch:
                    # Sample_name sera P\patT\time
                    sample_name = "P" + mtch.group("pat") + \
                        "T" + mtch.group("time")
                else:  # panic
                    raise Exception("Sample title do not follow the pattern.")
            s2gsm_file.write(sample_name+"\t"+values[sample_geo_acc_idx]+"\n")
            to_geo[sample_name] = values[sample_geo_acc_idx]
    s2gsm_file.close()

    # pegar arquivo GSE28405_non-normalized.txt.gz
    non_norm_fname = "data/geo_raw/supplemental_data/GSE28405/GSE28405_non-normalized.txt.gz"
    # descompacta-lo
    # ler e ignorar linhas ate chegar em ^ID_REF esse é o header
    detection_file = open(det_fname, "w")
    signal_file = open(signal_fname, "w")
    with gzip.open(non_norm_fname) as non_norm:
        header = ""
        for line in non_norm:
            line = line.decode("utf-8").strip()
            values = line.split("\t")
            if header:
                signal_values = [values[n] for n in non_det]
                signal_values.insert(0, values[0])
                signal_file.write("\t".join(signal_values)+"\n")
                det_values = [values[d] for d in det]
                det_values.insert(0, values[0])
                detection_file.write("\t".join(det_values)+"\n")
            if not header and values[0] == "ID_REF":
                header = values
                det = [i for i, v in enumerate(header) if v == 'Detection Pval']
                non_det = set(range(1, len(header))) - set(det)
                non_det = [n for n in non_det if header[n]]  # remove vazios

                signal_header = [to_geo[header[n]] for n in non_det]
                signal_header.insert(0, "ProbeName")
                signal_file.write("\t".join(signal_header)+"\n")
                det_header = [to_geo[header[d-1]] for d in det]
                det_header.insert(0, "ProbeName")
                detection_file.write("\t".join(det_header)+"\n")
    detection_file.close()
    signal_file.close()
    # descobrir o Sample_geo_accession de cada coluna
    # Imprimir tabelas de detection_pval e expression
