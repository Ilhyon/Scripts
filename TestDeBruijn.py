#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import networkx as nx
import matplotlib.pyplot as plt

def de_bruijn_ize(st, k):
    G = nx.Graph()
    for i in range( len(st) - k + 1 ):
        e = ( (st[i:i+k], st[i+1:i+k+1]) )
        G.add_node( st[i:i+k] )
        G.add_node( st[i+1:i+k+1] )
        if (st[i:i+k],st[i+1:i+k+1]) in list(G.edges):
            G[ st[i:i+k] ][ st[i+1:i+k+1] ]['weight'] += 1
        else:
            G.add_edge(st[i:i+k], st[i+1:i+k+1], weight=1 )
    # print(nx.get_edge_attributes(G,'weight'))
    # plt.subplot(121)
    # pos = nx.spring_layout(G, scale=2)
    # edge_labels = nx.get_edge_attributes(G,'weight')
    # nx.draw(G, with_labels=True)
    # nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels)
    # plt.show()
    pos = nx.spring_layout(G)
    edge_labels = nx.get_edge_attributes(G,'weight')
    nx.draw(G, pos, with_labels=True)
    nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels,arrows=True)
    plt.show()
    return G

def main():
    seq1 = 'GGGGGCCGGCUGGGAGGGCUGUCGGUGGGCCAGUCUGC'
    seq2 = 'GGGGCGGGGCCGGGAGGG'
    seq3 = 'TTGGGGCTGGCTGATGCTGGGGAGTGGGACTGTCTGCTGGGTTCAAGGTGGAAAGGGAAC'
    de_bruijn_ize(seq1, 3)
    de_bruijn_ize(seq2, 3)
    de_bruijn_ize(seq3, 3)

main()
# nx.test()
