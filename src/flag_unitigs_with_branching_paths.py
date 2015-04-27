import os, sys
import networkx as nx


gmlfn = sys.argv[1]

G= nx.read_gml(gmlfn)

exclude_utgs = set()
for node in G.nodes():
    if 	len(G.predecessors(node)) > 1 or len(G.successors(node)) > 1:
        exclude_utgs.add(G.node[node]['unitig'])

print "\n".join(list(exclude_utgs))
