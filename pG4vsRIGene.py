#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import re
import os
import argparse
import numpy as np
import pandas as pd
from pprint import pprint

def overlaps(interval1, interval2):
    """Computes the distance or overlap between two interval.

    Compute : min(ends) - max(starts). If > 0, return the number of bp of
    overlap, if 0, they are book-ended and if < 0 return the distance in
    bp between them.

    :param interval1: contain the start and the end of a OQs.
	:type interval1: list of int
    :param interval2: contain the start and the end of a gff feature.
	:type interval2: list of int

    :returns: min(ends) - max(starts)
    :rtype: int
    """
    return min(interval1[1], interval2[1]) - max(interval1[0], interval2[0])

def getLocationpG4(w, pG4, windowUp, windowDown):
    pG4coords = [ pG4['Start'], pG4['End'] ]
    # RIcoords = [ w['riExonStart_0base'], w['riExonEnd'] ]
    RIcoords = [ w -windowUp, w +windowDown]
    o = overlaps(pG4coords, RIcoords)
    if o > 0:
        return 1
    else:
        return 0

def readLineRI(w):
    line = {
    'newID' : w[0],
    'ID': w[2],
    'GeneID': w[3],
    'geneSymbol': w[4],
    'chr': w[5],
    'strand': w[6],
    'SplicingStart': int(w[7]),
    'SplicingEnd': int(w[8]),
    'upstreamES': int(w[9]),
    'upstreamEE': int(w[10]),
    'downstreamES': int(w[11]),
    'downstreamEE': int(w[12]),
    'PValue': w[20],
    'FDR': w[21],
    'significant': w[29]}
    return line

def read(filename, pG4, windowUp, windowDown, v, signi):
    with open(filename) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l and not l.startswith('newID2'):
                w = l.split('\t')
                w = readLineRI(w)
                if w['significant'] == signi:
                    if w['GeneID'] in pG4:
                        for G4 in pG4[ w['GeneID'] ]:
                            if w['strand'] == '+':
                                locUp = getLocationpG4(w['upstreamEE'], pG4[ w['GeneID'] ][G4], windowUp, windowDown)
                                if locUp != 0 :
                                    if w['geneSymbol'] in ['MAPK12', 'WDR90', 'SLC25A18', 'ULK3', 'IRF7']:
                                        print('upstream')
                                        print(v)
                                        print(w['geneSymbol'])
                                        print(w)
                                        print(pG4[ w['GeneID'] ][G4])
                                        print('----------')
                                locDown = getLocationpG4(w['downstreamES'], pG4[ w['GeneID'] ][G4], windowDown, windowUp)
                                if locDown != 0 :
                                    if w['geneSymbol'] in ['MAPK12', 'WDR90', 'SLC25A18', 'ULK3', 'IRF7']:
                                        print('downstream')
                                        print(v)
                                        print(w['geneSymbol'])
                                        print(w)
                                        print(pG4[ w['GeneID'] ][G4])
                                        print('----------')
                            else:
                                locDown = getLocationpG4(w['downstreamES'], pG4[ w['GeneID'] ][G4], windowUp, windowDown)
                                if locDown != 0 :
                                    if w['geneSymbol'] in ['MAPK12', 'WDR90', 'SLC25A18', 'ULK3', 'IRF7']:
                                        print('upstream')
                                        print(v)
                                        print(w['geneSymbol'])
                                        print(w)
                                        print(pG4[ w['GeneID'] ][G4])
                                        print('----------')
                                locUp = getLocationpG4(w['upstreamEE'], pG4[ w['GeneID'] ][G4], windowDown, windowUp)
                                if locUp != 0 :
                                    if w['geneSymbol'] in ['MAPK12', 'WDR90', 'SLC25A18', 'ULK3', 'IRF7']:
                                        print('downstream')
                                        print(v)
                                        print(w['geneSymbol'])
                                        print(w)
                                        print(pG4[ w['GeneID'] ][G4])
                                        print('----------')

def importpG4(filename):
    pG4 = {}
    with open(filename) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l and not l.startswith('Transcript') and not l.startswith('InfoG4ByTranscript'):
                w = l.split('\t')
                if w[9] not in pG4:
                    pG4[ w[9] ] = {}
                start = w[1].split(':')[1].split('-')[0]
                end = w[1].split(':')[1].split('-')[1]
                pG4[ w[9] ][ w[9]+'|'+w[1] ] = {'Start' : int(start),
                    'End' : int(end),
                    'Chr' : w[1].split(':')[0],
                    'Strand' : w[2],
                    'cGcC' : w[3],
                    'G4H' : w[4],
                    'G4nn' : w[6],
                    'Seq' : w[5],
                    'Bt' : w[8]}
    return pG4

def getpG4NearRI(path, windowUp, windowDown):
    pG4All = path + 'HS_pG4.csv'
    files = {'kunvRI' : path + 'RI/KUNV_RI.csv',
            'sinvRI' : path + 'RI/SINV_RI.csv',
            'zikvRI' : path + 'RI/ZIKV_RI.csv',
            'yvfRI' : path + 'RI/YFV_RI.csv'}
    pG4 = importpG4(pG4All)
    for v in files:
        # 1 = significant, 0 = non_significant
        read(files[v], pG4, windowUp, windowDown, v, '1')
        # read(files[v], pG4, windowUp, windowDown, v, '0')



def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-wUp', '--windowUp', default = 200)
    parser.add_argument ('-wDown', '--windowDown', default = 100)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    windowUp = arg.windowUp
    windowDown  = arg.windowDown
    getpG4NearRI(path, windowUp, windowDown)
