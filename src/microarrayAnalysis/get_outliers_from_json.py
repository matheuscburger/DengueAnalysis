#!/usr/bin/env python3
# -*- coding: utf8 -*-

"""Get Outliers from outliers.json inside aqm directory

Usage:
  get_outliers_from_json.py JSON --counts=<value>
  get_outliers_from_json.py (-h | --help)
  get_outliers_from_json.py --version

Options:
  -h --help                Show this screen.
  --version                Show version.
  JSON                     input json file
  --counts=<value>          Number of methods to consider a sample outlier
"""


__author__ = "Matheus Carvalho Bürger"
__email__ = "matheus.cburger@gmail.com"
__license__ = "GPL"


from docopt import docopt
import json


if __name__ == "__main__":
    args = docopt(__doc__, version='Get outliers from JSON')
    counts_min = int(args["--counts"])
    json_file = args["JSON"]
    js = json.load(open(json_file))
    outliers = [sample for sample, count in js['counts'].items()
                if int(count) >= counts_min]
    if len(outliers) > 0:
        print("\n".join(outliers))
