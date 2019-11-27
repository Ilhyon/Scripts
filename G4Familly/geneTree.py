#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import argparse
from Bio import Phylo

def pruneTree(tree, listToPrune):
    for clade in listToPrune:
        try:
            tree.prune(clade)
        except:
            tree = ''
    return tree

def getUnwantedGenes(tree, listSp):
    toPrune = []
    for clade in tree.get_terminals():
        # print(clade.taxonomy)
        taxonomyName = str(clade.taxonomy)
        tmp = [sp for sp in listSp if sp not in taxonomyName]
        tmp = list(set(tmp))
        if len(tmp) == len(listSp):
            toPrune.append(clade)
    return toPrune

def importListSp(filename):
    listSp = []
    with open(filename) as f:
        lines = f.read().splitlines()
        for l in lines:
            listSp.append(l)
    return listSp

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'geneTree')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

def main(listSp, directoryTrees, output):
    cpt = 0
    for path, dirs, files in os.walk(directoryTrees):
        # for each element of the directoryTrees to passed
        for filenameTree in files: # for each files
            inputTreefile = directoryTrees + '/' + filenameTree
            tree = Phylo.read(inputTreefile,'phyloxml')
            # print('++++++++++')
            # Phylo.draw_ascii(tree)
            listSpToPrune = getUnwantedGenes(tree, listSp)
            tree = pruneTree(tree, listSpToPrune)
            if tree and len(tree.get_terminals()) > 1:
                # Phylo.draw_ascii(tree)
                Phylo.write(tree, output + '/' + filenameTree, 'phyloxml')
            else:
                cpt += 1
    print(cpt)

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    listSp = importListSp(path + 'Genomes/Allsp/listSp.txt')
    directoryTrees = path + 'Tree/EnsemblTree'
    output = path + 'Tree/PrunedTree'
    main(listSp, directoryTrees, output)
