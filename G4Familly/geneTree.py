#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

from Bio import Phylo

treeFile = '/home/anais/Documents/Data/Tree/Compara.97.protein_default.tree.phyloxml.xml.1000566-1006055/protein_default.tree.1000566.phyloxml.xml'
tree = Phylo.read(treeFile,'phyloxml')
Phylo.draw_ascii(tree)
print(tree)
