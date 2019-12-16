#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import os
import argparse
from Bio import SeqIO
from pprint import pprint
import rF

def reverseSequence(Sequence) :
    """ Reverse complement a DNA sequence.

    :param Sequence: DNA sequence that will be reversed.
    :type Sequence: string

    :returns: Sequence, the initial DNA sequence but reverse complemented.
    :rtype: string
    """
    reverse = ""
    nucleotides = {'A' : 'T',
                    'T' : 'A',
                    'C' : 'G',
                    'G' : 'C'}
    for n in Sequence:
        if n in nucleotides:
            tmp = nucleotides[n]
        else :
            tmp = n # in some sequences there is many N or other letter
    reverse += tmp
    Sequence = reverse[::-1]
    return Sequence

def getFasta(path, dicoGeneTree, dicoGeneList):
    dicoFastaTree = {}
    for tree in dicoGeneTree:
        dicoFastaTree[tree] = {}
        for gene in dicoGeneTree[tree]:
            spId = dicoGeneList[gene].split('-')[0]
            sp = dicoGeneList[gene].split('-')[1]
            filename = path + 'Genomes/' + sp + '/' + spId + '_gene_unspliced.txt'
            geneUnspliced = SeqIO.to_dict(SeqIO.parse(filename, 'fasta'))
            print(record_dict[gene])
            dicoFastaTree[tree][gene] = record_dict[gene]


def importGeneList(geniListDir):
    dicoGeneListbySp = {}
    dicoSpId = rF.getDicoSpId()
    for path, dirs, files in os.walk(geniListDir):
        # for each element of the directory to passed
        for filename in files: # for each files
            inputfile = geniListDir + '/' + filename
            spId = filename.split('_')[0]
            if spId in dicoSpId:
                sp = dicoSpId[spId]
                with open(inputfile) as f:
                    content = f.read()
                    lines = content.split('\n')
                    for l in lines:
                        dicoGeneListbySp[l] = spId + '-' + sp
    return dicoGeneListbySp

def importGeneTree(treeGeneListDir):
    dicoGeneTree = {}
    for path, dirs, files in os.walk(treeGeneListDir):
        # for each element of the directory to passed
        for filename in files: # for each files
            inputfile = treeGeneListDir + '/' + filename
            treeId = filename.split('_')[1]
            dicoGeneTree[treeId] = []
            with open(inputfile) as f:
                content = f.read()
                lines = content.split('\n')
                for l in lines:
                    dicoGeneTree[treeId].append(l)
    return dicoGeneTree

def main(path):
    treeGeneListDir = path + 'Tree/PrunedGeneList/'
    geniListDir = path + 'Homology/'
    outputDir = path + 'Tree/Fasta/'
    dicoGeneTree = importGeneTree(treeGeneListDir)
    dicoGeneList = importGeneList(geniListDir)
    dicoFasta = getFasta(path, dicoGeneTree, dicoGeneList)

def build_arg_parser():
    GITDIR = os.getcwd()+'/'
    parser = argparse.ArgumentParser(description = 'getFastaTree')
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    main(path)
