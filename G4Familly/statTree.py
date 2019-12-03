#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import argparse
from Bio import Phylo
from pprint import pprint

def noRedondance(tree, dicoSp, filename):
    for clade in tree.get_terminals():
        taxonomyName = str(clade.taxonomy)
        if taxonomyName not in dicoSp:
            dicoSp[taxonomyName] = []
        dicoSp[taxonomyName].append(filename)
    return dicoSp

def Allsp(tree, filename):
    listSp = []
    for clade in tree.get_terminals():
        taxonomyName = str(clade.taxonomy)
        listSp.append(taxonomyName)
    if len(set(listSp)) >= 60:
        print(filename)
        Phylo.draw_ascii(tree)
        print('-----------------------------------------------------------------')

def Redondance(tree, dicoSp):
    for clade in tree.get_terminals():
        taxonomyName = str(clade.taxonomy)
        if taxonomyName not in dicoSp:
            dicoSp[taxonomyName] = 0
        dicoSp[taxonomyName] += 1
    return dicoSp

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'pruneTree')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

def main(directoryTrees):
    dicoSpRedondance = {}
    dicoSpNoRedondance = {}
    for path, dirs, files in os.walk(directoryTrees):
        # for each element of the directoryTrees to passed
        for filenameTree in files: # for each files
            inputTreefile = directoryTrees + '/' + filenameTree
            tree = Phylo.read(inputTreefile,'phyloxml')
            Allsp(tree, filenameTree)
            dicoSpRedondance = Redondance(tree, dicoSpRedondance)
            dicoSpNoRedondance = noRedondance(tree, dicoSpNoRedondance, filenameTree)
    for sp in dicoSpNoRedondance:
        dicoSpNoRedondance[sp] = len(set(dicoSpNoRedondance[sp]))
    pprint(dicoSpRedondance)
    pprint(dicoSpNoRedondance)

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    directoryTrees = path + 'Tree/PrunedTree'
    main(directoryTrees)
