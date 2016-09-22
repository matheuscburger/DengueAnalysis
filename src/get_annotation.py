#!/usr/bin/env python3
# vim:fileencoding=utf8

"""Get Sample Annotation

Usage:
  get_annotation.py --reannotation=<file> --platform=<GPL>
  get_annotation.py (-h | --help)
  get_annotation.py --version

Options:
  -h --help                Show this screen.
  --version                Show version.
  --reannotation=<file>    Probe Annotation file
  --platform=<GPL>         Platform ID

"""

__author__ = "Matheus Carvalho Bürger"
__email__ = "matheus.cburger@gmail.com"
__license__ = "GPL"

from docopt import docopt

if __name__ == "__main__":
    args = docopt(__doc__, version='Get annotation')
    reannotation_filename = args["--reannotation"]
    platform = args["--platform"]
    with open(reannotation_filename) as reannotation_file:
        header = reannotation_file.readline().strip("\n").split("\t")
        nannot_idx = header.index("NumAnnot")
        hits_idx = header.index("Hits")
        platform_idx = header.index("Platform")
        gene_idx = header.index("Gene")
        probe_idx = header.index("Probe")
        # Header
        print("ProbeName\tSymbol")
        for line in reannotation_file:
            line = line.strip("\n")
            values = [v.strip() for v in line.split("\t")]
            if int(values[hits_idx]) == 1 and int(values[nannot_idx]) == 1 and \
                    values[platform_idx] == platform:
                print("%s\t%s" % (values[probe_idx], values[gene_idx]))
