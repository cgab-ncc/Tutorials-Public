# Creating an undirected graph
import networkx as nx
G1 = nx.Graph()

# Add nodes
G1.add_node(1)

# Add multiple nodes
G1.add_nodes_from([2,3,4])

# Add edges
G1.add_edge(1,2)
G1.add_edge(1,3)

# Add multiple edges
G1.add_edges_from([(2,3), (3,4)]) 

# Stats
print("G1 has "  + str(G1.number_of_nodes()) + " nodes")
print("G1 has " + str(G1.number_of_edges()) + " edges")

# See nodes node 1 is connected to
print("Nodes node 1 is connected to : " + str(G1[1]))

# Visualize graph
import matplotlib.pyplot as plt

nx.draw(G1, with_labels=True, node_color='y')
# nx.draw_random(G1)
# nx.draw_circular(G1)
# nx.draw_spectral(G1)
plt.show()

# Add attributes to nodes
G2 = nx.Graph()
G2.add_node(1, gene_exp='5.17')
G2.add_node(2, gene_exp='0.53')
G2.add_node(3, gene_exp='1.81')

# Add weight to edge
G2.add_edge(1, 2, weight=4.7)
print("Weight of edge(1,2) : " + str(G2[1][2]))

# Update edge weight
G2[1][2]['weight'] = 2.9
print("Updated weight of edge(1,2) : " + str(G2[1][2]))

# Creating a directed graph
DG=nx.DiGraph()
DG.add_weighted_edges_from([(1,2,0.5), (3,1,0.75), (1,4,0.33)]) # (from, to, weight)
DG.out_degree(1,weight='weight')
DG.degree(1,weight='weight')
print("Successors of node 1 : " + str(DG.successors(1)))
print("Neighbors of node 1 : " + str(DG.neighbors(1)))
print("In degree of node 1 : " +str(DG.in_degree(1)))
print("Out degree of node 1 : " +str(DG.out_degree(1)))

# Visualizing a directed graph
import pylab
edge_labels=dict([((u,v,),d['weight'])for u,v,d in DG.edges(data=True)])
pos=nx.spring_layout(DG)
nx.draw_networkx_edge_labels(DG,pos,edge_labels=edge_labels)
nx.draw(DG,pos,edge_cmap=plt.cm.Reds, with_labels=True, node_color='y')
pylab.show()