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

def getUnwantedGenes(tree, listSp, listMissing):
    toPrune = []
    listGene = []
    for clade in tree.get_terminals():
        taxonomyName = str(clade.taxonomy)
        #the follwing is just to check the real 'name' of a species if it is
        #missing in my results
        [print(taxonomyName) for sp in listMissing if sp in taxonomyName]
        tmp = [sp for sp in listSp if sp not in taxonomyName]
        tmp = list(set(tmp))
        if len(tmp) == len(listSp):
            toPrune.append(clade)
        else:
            listGene.append(str(clade.name))
    listGene = list(set(listGene))
    return toPrune, listGene

def importListSp(filename):
    listSp = []
    with open(filename) as f:
        lines = f.read().splitlines()
        for l in lines:
            listSp.append(l)
    return listSp

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'pruneTree')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

def main(listSp, directoryTrees, outputTree, outputGeneList):
    cpt = 0
    listMissing = []
    for path, dirs, files in os.walk(directoryTrees):
        # for each element of the directoryTrees to passed
        for filenameTree in files: # for each files
            inputTreefile = directoryTrees + '/' + filenameTree
            tree = Phylo.read(inputTreefile,'phyloxml')
            listSpToPrune, geneListPruned = getUnwantedGenes(tree, listSp, listMissing)
            tree = pruneTree(tree, listSpToPrune)
            if tree and len(tree.get_terminals()) > 2:
                Phylo.write(tree, outputTree + '/' + filenameTree, 'phyloxml')
                listName = filenameTree.split('.')[2]
                output = open(outputGeneList + '/geneList_' + listName ,'w')
                output.write('\n'.join(geneListPruned))
            else:
                cpt += 1
    print(cpt)

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    listSp = importListSp(path + 'Genomes/Allsp/listSp.txt')
    directoryTrees = path + 'Tree/EnsemblTree'
    outputTree = path + 'Tree/PrunedTree'
    outputGeneList = path + 'Tree/PrunedGeneList'
    main(listSp, directoryTrees, outputTree, outputGeneList)
