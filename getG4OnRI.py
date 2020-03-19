#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import re
import os
import argparse
import pandas as pd
from pprint import pprint

def getListTrIntron(filename):
    df = pd.read_csv(filename, sep='\t')
    df = df[ df.significant == 1]
    geneId = list(df.GeneID)
    return geneId

def getNbG4PerVirus(path):
    dfpG4v = {}
    pG4All = path + 'HS_pG4.csv'
    viruses = {'kunvRI' : path + '200319_KUNV_RI.csv',
        'sinvRI' : path + '200319_SINV_RI.csv',
        'zikvRI' : path + '200319_ZIKV_RI.csv',
        'yvfRI' : path + '200319_YFV_RI.csv'}
    dfpG4 = pd.read_csv(pG4All, sep='\t')
    for v in viruses:
        geneId = getListTrIntron(viruses[v])
        dfpG4v[v] = dfpG4[dfpG4['Gene'].isin(geneId)]
        del dfpG4v[v]['Transcript']
        dfpG4v[v] = dfpG4v[v].drop_duplicates(subset=None, keep='first', inplace=False)
        # print(v)
        # print(len(dfpG4v[v]))

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

def getLocationpG4(w, pG4, window):
    pG4coords = [ pG4['Start'], pG4['End'] ]
    RIcoords = [ w['riExonStart_0base'], w['riExonEnd'] ]
    # RIcoords = [ w['upstreamEE'], w['downstreamES'] ]
    o = overlaps(pG4coords, RIcoords)
    if o > 0:
        if o == (pG4coords[1] - pG4coords[0]):
            location = 'pG4InRI'
        elif pG4coords[0] < RIcoords[0]:
            location = "pG4overlap5RI'"
        else:
            location = "pG4overlap3RI'"
    elif o == 0:
        location = 'pG4TouchRI'
    elif o < 0 and o > - window:
        location = 'pG4nearRI'
    else:
        location = ''
    return location


def readLineRI(w):
    line = {'newID2': w[0],
    'newID': w[1],
    'ID': w[2],
    'GeneID': w[3],
    'geneSymbol': w[4],
    'chr': w[5],
    'strand': w[6],
    'riExonStart_0base': int(w[7]),
    'riExonEnd': int(w[8]),
    'upstreamES': int(w[9]),
    'upstreamEE': int(w[10]),
    'downstreamES': int(w[11]),
    'downstreamEE': int(w[12]),
    'ID_1': w[13],
    'IC_SAMPLE_1': w[14],
    'SC_SAMPLE_1': w[15],
    'IC_SAMPLE_2': w[16],
    'SC_SAMPLE_2': w[17],
    'IncFormLen': w[18],
    'SkipFormLen': w[19],
    'PValue': w[20],
    'FDR': w[21],
    'IncLevel1': w[22],
    'misc1': w[23],
    'misc2': w[24],
    'IncLevel2': w[25],
    'misc3': w[26],
    'misc4': w[27],
    'IncLevelDifference': w[28],
    'significant': w[29]}
    return line

def readRI(filename, pG4, window):
    RIpG4 = {}
    cpt = 0
    with open(filename) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l and not l.startswith('newID2'):
                w = l.split('\t')
                w = readLineRI(w)
                if w['GeneID'] in pG4:
                    for G4 in pG4[ w['GeneID'] ]:
                        loc = getLocationpG4(w, pG4[ w['GeneID'] ][G4], window)
                        if loc:
                            RIpG4[cpt] = {}
                            RIpG4[cpt].update(w)
                            RIpG4[cpt].update( pG4[ w['GeneID'] ][G4])
                            RIpG4[cpt]['Location'] = loc
                            cpt += 1
    pprint(RIpG4)

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
                if w[2] == '1':
                    start = w[1].split(':')[1].split('-')[0]
                    end = w[1].split(':')[1].split('-')[1]
                else:
                    end = w[1].split(':')[1].split('-')[0]
                    start = w[1].split(':')[1].split('-')[1]
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

def getpG4NearRI(path, window):
    pG4All = path + 'HS_pG4.csv'
    viruses = {'kunvRI' : path + '200319_KUNV_RI.csv',
        'sinvRI' : path + '200319_SINV_RI.csv',
        'zikvRI' : path + '200319_ZIKV_RI.csv',
        'yvfRI' : path + '200319_YFV_RI.csv'}
    pG4 = importpG4(pG4All)
    readRI(path + '200319_KUNV_RI.csv', pG4, window)


def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-w', '--window', default = 100)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    window = arg.window
    getpG4NearRI(path, window)
