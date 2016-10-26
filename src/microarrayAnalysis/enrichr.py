#!/usr/bin/env python3
# vim:fileencoding=utf8

"""Remove columns from a table (TSV)

Usage:
  enrichr.py --input=<file> --output=<file> (--gs=<geneset>...)
  enrichr.py (-h | --help)
  enrichr.py --version

Options:
  -h --help                 Show this screen.
  --version                 Show version.
  --output=<file>           output file
  --input=<file>            input file, each gene in one row
"""

__author__ = "Matheus Carvalho Bürger"
__email__ = "matheus.cburger@gmail.com"
__license__ = "GPL"

import json
import requests
from docopt import docopt

def getResultsFile(output, genesets, userlistid,
                   enrichr_url="http://amp.pharm.mssm.edu/",
                   query_string="Enrichr/export?userListId=%s&filename=%s&backgroundType=%s"):
    flagHeader = False
    with open(output, "w") as f:
        for gs in genesets:
            url = enrichr_url + query_string % \
                (userlistid, output, gs)
            response = requests.get(url, stream=True)
            lines = response.text.strip("\n").split("\n")
            header = "GeneSet\t"+lines.pop(0)
            if not flagHeader:
                print(header, file=f)
                flagHeader = True
            lines = [gs+"\t"+l+"\n" for l in lines]
            f.writelines(lines)


def addList(genes, description, enrichr_url="http://amp.pharm.mssm.edu/",
            addList_url="Enrichr/addList"):
    payload = {
        "list": (None, genes_str),
        "description": (None, description)
    }
    response = requests.post(enrichr_url + addList_url, files=payload)
    if not response.ok:
        raise Exception("Error analyzing gene list")
    return(json.loads(response.text))


if __name__ == "__main__":
    args = docopt(__doc__, version='1.0')

    genes_fname = args["--input"]
    outfname = args["--output"]
    genesets = args["--gs"]

    with open(genes_fname) as genes_fh:
        genes_str = genes_fh.read()

    list_dict = addList(genes_str, "")
    getResultsFile(outfname, genesets, list_dict["userListId"])
