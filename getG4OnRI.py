#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import re
import os
import argparse
import pandas as pd
from pprint import pprint

def getGeneTr(filename):
    dico = {}
    with open(filename) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l:
                w = l.split('|')
                dico[ w[0] ] = w[1]
    return dico

def readpG4(filename, genTrdico):
    linesN = []
    with open(filename) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l:
                w = l.split('\t')
                if w[0] in genTrdico:
                    gene = genTrdico[ w[0] ]
                else:
                    gene = 'Gene'
                w.append(gene)
                linesN.append( '\t'.join(w) )
        output = open('HS_pG4.csv', "w")
        for l in linesN:
            output.write(l + "\n")

def getListTrIntron(filename):
    df = pd.read_csv(filename, sep='\t')
    geneId = list(df.GeneID)
    return geneId

def main(path):
    dfpG4v = {}
    pG4All = path + 'HS_All_G4InTranscript.txt'
    genTr = path + 'HS_transcript_unspliced_All_Index.txt'
    viruses = {'kunvRI' : path + '200319_KUNV_RI.csv',
        'sinvRI' : path + '200319_SINV_RI.csv',
        'zikvRI' : path + '200319_ZIKV_RI.csv',
        'yvfRI' : path + '200319_YFV_RI.csv'}
    genTrdico = getGeneTr(genTr)
    readpG4(pG4All, genTrdico)
    # dfpG4 = pd.read_csv(pG4All, sep='\t')
    # dfpG4['Gene'] = dfpG4.Transcript.apply(addGene)
    # dfpG4.to_csv(path_or_buf=path+'HS_pG4.csv', header=True, index=None, sep='\t')
    # for v in viruses:
    #     geneId = getListTrIntron(viruses[v])
    #     dfpG4v[v] = dfpG4[dfpG4['Gene'].isin(geneId)]
    #     print(v)
    #     print(len(dfpG4v[v]))




def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    main(path)
