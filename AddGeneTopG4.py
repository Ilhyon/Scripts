#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os
from pprint import pprint

if __name__ == '__main__':
    GITDIR = os.getcwd()+'/'
    pG4All = GITDIR + 'NewHSpG4.csv'
    index = '/home/anais/Documents/G4HumanTranscriptome/Data/transcriptType/transcriptType_All.txt'
    d = {}
    pG4 = []
    with open(index) as f:
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l and not l.startswith('Transcript') and not l.startswith('InfoG4ByTranscript'):
                w = l.split('\t')
                d[w[1]] = w[0]
    with open(pG4All) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l and not l.startswith('Transcript') and not l.startswith('InfoG4ByTranscript'):
                w = l.split('\t')
                if len(w) == 9:
                    w.append(d[w[0]])
                w = '\t'.join(w)
                pG4.append(w)
    outputF = open(GITDIR+'HS_pG4.csv', "w")
    outputF.write( '\n'.join(pG4) )
    outputF.close()
