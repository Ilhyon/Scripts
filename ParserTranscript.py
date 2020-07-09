#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os
import argparse
import numpy as np
import pandas as pd
from pprint import pprint

def readTr(filename):
    dico = {}
    with open(filename) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l :
                w = l.split('|')
                if w[3] == '1':
                    w[3] = '+'
                else:
                    w[3] = '-'
                chrStrand = w[2]+'|'+w[3]
                if chrStrand not in dico:
                    dico[chrStrand] = {}
                exon = w[5].split(';')
                for e in exon:
                    if e not in dico[chrStrand]:
                        dico[chrStrand][e] = []
                    dico[chrStrand][e].append(w[0])
    return dico

def main(path):
    trAll = path + 'HS_transcript_unspliced_All.txt'
    files = ['kunv', 'sinv', 'zikv', 'yvf']
    dicoAllTr = readTr(trAll)
    for v in files:
        newF = []
        with open(path+v+'_RI1New.csv') as f: # file opening
            content = f.read()
            lines = content.split('\n')
            for l in lines:
                tr1 = []
                tr2 = []
                w = l.split('\t')
                if w[2] == '-':
                    E1E = str(int(w[9])+1)
                    E1S = str(int(w[10]))
                    E2E = str(int(w[11])+1)
                    E2S = str(int(w[12]))
                    chrStrand = w[3]+'|'+w[2]
                    if E1S+'-'+E1E in dicoAllTr[chrStrand]:
                        tr1 = dicoAllTr[chrStrand][ E1S+'-'+E1E ]
                    else:
                        print('tr1')
                        print(E1S+'-'+E1E)
                    if E2S+'-'+E2E in dicoAllTr[chrStrand]:
                        tr2 = dicoAllTr[chrStrand][ E2S+'-'+E2E ]
                    else:
                        print('tr2')
                        print(E2S+'-'+E2E)
                    if tr1 and tr2:
                        commonTr = list(set(tr1).intersection(tr2))
                    else:
                        commonTr = []
                    w.extend(commonTr)
                    w = '\t'.join(w)
                    newF.append(w)
                else:
                    E1S = str(int(w[9])+1)
                    E1E = str(int(w[10]))
                    E2S = str(int(w[11])+1)
                    E2E = str(int(w[12]))
                    chrStrand = w[3]+'|'+w[2]
                    if E1S+'-'+E1E in dicoAllTr[chrStrand]:
                        tr1 = dicoAllTr[chrStrand][ E1S+'-'+E1E ]
                    else:
                        print('tr1')
                        print(E1S+'-'+E1E)
                    if E2S+'-'+E2E in dicoAllTr[chrStrand]:
                        tr2 = dicoAllTr[chrStrand][ E2S+'-'+E2E ]
                    else:
                        print('tr2')
                        print(E2S+'-'+E2E)
                    if tr1 and tr2:
                        commonTr = list(set(tr1).intersection(tr2))
                    else:
                        commonTr = []
                    w.extend(commonTr)
                    w = '\t'.join(w)
                    newF.append(w)
        outputF = open(path+v+'_RI1TESTtranscript.csv', "w")
        outputF.write( 'Location\tGeneSymbol\tStrand\tchr\tStartEvent\tEndEvent\tStartpG4\tEndpG4\tpG4Sequence\tE1S\tE1E\tE2S\tE2E\tTr\n' )
        outputF.write( '\n'.join(newF) )
        outputF.close()


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
