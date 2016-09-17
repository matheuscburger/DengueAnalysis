#!/usr/bin/env python3
# -*- coding: utf8 -*-

"""Remove columns from a table (TSV)

Usage:
  remove_columns.py INPUT [--columns=<value>...]
  remove_columns.py (-h | --help)
  remove_columns.py --version

Options:
  -h --help                 Show this screen.
  --version                 Show version.
  INPUT                     input file
  --counts=<value>          Number of methods to consider a sample outlier
"""


__author__ = "Matheus Carvalho Bürger"
__email__ = "matheus.cburger@gmail.com"
__license__ = "GPL"


from docopt import docopt


def remove_idx(vals, idx):
    return([v for i, v in enumerate(vals) if i not in idx])

if __name__ == "__main__":
    args = docopt(__doc__, version='Get outliers from JSON')
    columns = args["--columns"]
    filename = args["INPUT"]
    with open(filename) as fh:
        header_line = fh.readline().strip("\n")
        header = header_line.split("\t")
        idx = [i for i, v in enumerate(header) if v in columns]
        print("\t".join(remove_idx(header, idx)))
        for line in fh:
            line = line.strip("\n")
            values = line.split("\t")
            print("\t".join(remove_idx(values, idx)))
