#!/usr/bin/env python3
# vim:fileencoding=utf8

__author__ = "Matheus Carvalho Bürger"
__email__ = "matheus.cburger@gmail.com"
__license__ = "GPL"

import graph_tool.all
import sys

if __name__ == "__main__":
    graphml = sys.argv[1]
    prefix = sys.argv[2]
    g = graph_tool.all.load_graph(graphml)
    clusters = [k for k in g.vertex_properties.keys() if not k.startswith("_")
                and k != "name"]
    pos = graph_tool.all.sfdp_layout(g)

    bool_filt = []
    view_filt = []
    for i in range(3):
        bool_filt.append([e >= i+1 for e in g.edge_properties["Sum"]])
        view_filt.append(graph_tool.all.GraphView(g, efilt=bool_filt[i]))

    for c in clusters:
        for i in range(3):
            graph_tool.all.graph_draw(view_filt[i], pos=pos,
                                      bg_color=[1, 1, 1, 1],
                                      output=prefix + c + "_get" + str(i+1) +
                                      ".png",
                                      vertex_fill_color=g.vertex_properties[c],
                                      edge_pen_width=g.edge_properties["Sum"])
