#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import os
import argparse
from Bio import SeqIO
from pprint import pprint
import rF

def writeFasta(outDir, dicoFasta, tree):
    """From a dictionary {id : seq}, write a fasta file.

    :param outDir: name of the directory where the output file need to be writen.
    :type outDir: string
    :param dicoFasta: {id : seq}
    :type dicoFasta: dictionary
    :param chr: name of the chromosome.
    :type chr: string
    """
    output = open(outDir +'FastaTree_' + tree + '.fas', "w")
    for id in dicoFasta:
        output.write('>' + id + "\n")
        nbLine = math.ceil( float( len(dicoFasta[id]) ) / 60 )
        cpt1 = 0
        cpt2 = 60
        for i in range(0,int(nbLine)) :
            seq = str(dicoFasta[id])
            output.write(str(seq[cpt1:cpt2]) + "\n")
            # to have a new line after 60 characters
            cpt1 += 60
            cpt2 += 60
    output.close()

def getFasta(path, dicoGeneTree, dicoGeneList, outDir):
    for tree in dicoGeneTree:
        dicoFastaTmp = {}
        for gene in dicoGeneTree[tree]:
            spId = dicoGeneList[gene].split('-')[0]
            sp = dicoGeneList[gene].split('-')[1]
            filename = path + 'Genomes/' + sp + '/' + spId + '_gene_unspliced.txt'
            geneUnspliced = SeqIO.to_dict(SeqIO.parse(filename, 'fasta'))
            parsedGeneUnspliced = {}
            infoGene = {}
            for seqId in geneUnspliced:
                parsedGeneUnspliced[ seqId.split(' ')[0] ] = geneUnspliced[seqId]
                infoGene[ seqId.split(' ')[0] ] = seqId
            dicoFastaTmp[ parsedGeneUnspliced[gene].description ] = parsedGeneUnspliced[gene].seq
        writeFasta(outDir, dicoFasta, tree)

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
    getFasta(path, dicoGeneTree, dicoGeneList, outputDir)

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
