#!/usr/bin/env python
# -*- coding: utf8 -*-

"""Get Sample Annotation

Usage:
  get_sample_annot.py [--out=<out-file>] <geo-id> ...
  get_sample_annot.py (-h | --help)
  get_sample_annot.py --version

Options:
  -h --help                Show this screen.
  --version                Show version.
  -o=<file> --out=<file>   output file

"""


__author__ = "Matheus Carvalho Bürger"
__email__ = "matheus.cburger@gmail.com"
__license__ = "GPL"


import re
import sys
import urllib.request
from Bio import Entrez
Entrez.email = "sysbio-usp@googlegroups.com"
from docopt import docopt


def create_dict(mapping):
    d = dict()
    for k, v in mapping:
        if k in d.keys():
            d[k] += "; "
        else:
            d[k] = ""
        d[k] += v
    return d


def get_accessions(accessions):
    retmax = int(1e9)   # parameter to retrieve all results
    acc_re = re.compile("Accession:\s*(?P<acc>\S+)\s+ID:\s*(?P<id>\d+)")
    accs = []
    search_acc = "("+" OR ".join(['"'+g+'"[Accession]' for g in accessions])+")"
    search = search_acc + 'AND "gsm"[Filter]'
    handle = Entrez.esearch(db="gds", term=search, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    if int(record["Count"]) > 0 and int(record['RetMax']) < retmax:
        idlist = [str(rid) for rid in record["IdList"]]
        handle = Entrez.efetch(db="gds", id=idlist)
        record = handle.read()
        for acc_mtch in acc_re.finditer(record):
            accs.append(acc_mtch.group("acc"))
    else:
        pass  # error: no results or too much results
    return(accs)


def get_sample_info_old(sample):
    block_init_re = re.compile("^\^SAMPLE\s*=\s*(?P<sample>\S+)")
    soft_re = re.compile("^\!(\S+)\s*=\s*(.*)$", flags=re.MULTILINE)
    url_base = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
    url = url_base+"acc=%s&targ=self&form=text&view=brief"
    response = urllib.request.urlopen(url % sample)
    soft = response.read().decode().replace("\r", "")
    match_init = block_init_re.match(soft)
    if match_init:
        if match_init.group("sample") == sample:
            soft_match = soft_re.findall(soft)
            sample_info = create_dict(soft_match)
        else:
            pass  # error: wrong SOFT file
    else:
        pass  # error: file should be in SOFT format
    return sample_info


def get_soft(study):
    url_base = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
    url = url_base+"acc=%s&targ=gsm&form=text&view=brief"
    response = urllib.request.urlopen(url % study)
    soft = response.read().decode().replace("\r", "")
    return(soft)


def get_sample_info(soft):
    block_re = re.compile(r"""^\^SAMPLE\s*=\s* # block init
                          (?P<sample>\S+)      # sample name
                          \n
                          (?P<info>(\!.+\n)+)  # sample information
                          """, flags=re.VERBOSE+re.MULTILINE)
    field_re = re.compile("^\!(\S+)\s*=\s*(.*)$", flags=re.MULTILINE)
    field_iter = block_re.finditer(soft)
    for block in field_iter:
        # sample = block.group("sample")
        info = block.group("info")
        fields_match = field_re.findall(info)
        sample_info = create_dict(fields_match)
        yield sample_info


if __name__ == "__main__":
    args = docopt(__doc__, version='Get sample annot 0.1')
    if args["--out"]:
        out = open(args["--out"], "w")
    else:
        out = sys.stdout
    studies = args["<geo-id>"]
    samples = get_accessions(studies)
    int_cols = ["Sample_series_id", "Sample_geo_accession",
                "Sample_platform_id", "Sample_supplementary_file",
                "Sample_title", "Sample_source_name_ch1",
                "Sample_characteristics_ch1"]
    print("\t".join(int_cols), file=out)
    for study in studies:
        soft = get_soft(study)
        for sample_info in get_sample_info(soft):
            line = [sample_info.get(cols, "NA") for cols in int_cols]
            print("\t".join(line), file=out)
