#!/usr/bin/env python3
# vim:fileencoding=utf8

"""Get modules (connected components) in a undirected graph represented by
adjacency list in a table (TSV file).

Usage:
  get_modules.py --input=<file> --from-col=<value> --to-col=<value> [--output=<file>] [(--filter-col=<value> --filter-val=<value>)]
  get_modules.py (-h | --help)
  get_modules.py --version

Options:
  -h --help                 Show this screen.
  --version                 Show version.
  --output=<file>           output file
  --input=<file>            input file
  --from-col=<value>        column containing a node label (from)
  --to-col=<value>          column containing a node label (to)
  --filter-col=<value>      numeric coulumn to filter relations
  --filter-val=<value>      numeric value, if filter column is below this value that row will not be considered
"""

__author__ = "Matheus Carvalho Bürger"
__email__ = "matheus.cburger@gmail.com"
__license__ = "GPL"

from docopt import docopt
import sys

class Graph():
    """
    Class that represents a graph
    """
    def __init__(self):
        """
        Initializer method
        It creates an attribute nodes that is a dictionary containing all nodes
        """
        self.nodes = {}

    def add_edge(self, labelA, labelB):
        """
        Adds edges to the graph, based on the node's name
        If the node does not exists it is created
        """
        nodeA = self.get_node(labelA)
        nodeB = self.get_node(labelB)
        nodeA.add_adjacent(nodeB)
        nodeB.add_adjacent(nodeA)

    def get_node(self, label):
        """
        Method that receives a node's name and returns the correspondent node
        If the node does not exists, it is created
        """
        if label not in self.nodes.keys():
            node = Node(label)
            self.nodes[label] = node
        else:
            node = self.nodes[label]
        return(node)


class Node():
    """
    Class that represents a node
    """
    def __init__(self, label):
        """
        Initializer method of the node
        It creates an attribute label and
        adjacents that is a dicionary containing all neighbors
        """
        self.label = label
        self.adjacents = {}

    def add_adjacent(self, node):
        """
        Adds a node to the adjacents dictionary
        """
        self.adjacents[node.label] = node

    def __str__(self):
        return(list(self.adjacents.keys()).__str__())

    def __repr__(self):
        return(self.__str__())


def depth_search_rec(node, mod, to_visit, visited, g):
    """
    Depth first search (recursive)
    It does not works with huge graphs
    """
    visited.add(node.label)
    mod.append(node.label)
    to_visit.extend([i for i in node.adjacents.keys() if i not in visited and
                     i not in to_visit])
    if len(to_visit) != 0:
        next_node = g.nodes[to_visit.pop()]
        depth_search(next_node, mod, to_visit, visited, g)

def depth_search(first, visited, g):
    """
    Busca em profundidade iterativa
    Depth first search
    Returns a list representing a module with all names of connected nodes
    """
    mod = []
    to_visit = [first]
    while len(to_visit) != 0:
        node = g.nodes[to_visit.pop()]
        visited.add(node.label)
        mod.append(node.label)
        to_visit.extend([i for i in node.adjacents.keys() if i not in visited and
                        i not in to_visit])
    return(mod)

def get_modules(g):
    """
    Receives a graph and returns modules i.e. connected components of a graph
    """
    modules = []
    visited = set()
    all_labels = list(g.nodes.keys())
    to_visit = []
    for current in all_labels:
        if current in visited:
            continue
        mod = depth_search(current, visited, g)
        modules.append(mod)
    return(modules)


if __name__ == "__main__":
    args = docopt(__doc__, version='1.0')

    in_fname = args["--input"]
    if args["--output"]:
        out_fname = args["--output"]
        out = open(out_fname, "w")
    else:
        out = sys.stdout
    from_col = args["--from-col"]
    to_col = args["--to-col"]
    filter_col = args["--filter-col"]
    if filter_col:
        filter_val = float(args["--filter-val"])

    graph = Graph()
    with open(in_fname) as cem:
        header = cem.readline().strip("\n").split("\t")
        g1_idx = header.index(from_col)
        g2_idx = header.index(to_col)
        if filter_col:
            filter_idx = header.index(filter_col)
        for i, line in enumerate(cem):
            line = line.strip("\n")
            values = line.split("\t")
            if filter_col and float(values[filter_idx]) < filter_val:
                continue
            graph.add_edge(values[g1_idx], values[g2_idx])
    modules = get_modules(graph)
    modules.sort(key = lambda x: len(x), reverse=True)
    for mod in modules:
        print(*sorted(mod), sep="\t", file=out)
    if args["--output"]:
        out.close()
