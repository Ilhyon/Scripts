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

def getLocationpG4(w, pG4, window):
    pG4coords = [ pG4['Start'], pG4['End'] ]
    # RIcoords = [ w['riExonStart_0base'], w['riExonEnd'] ]
    RIcoords = [ w['upstreamEE'] -100, w['downstreamES'] +100]
    o = overlaps(pG4coords, RIcoords)
    if o > 0:
        return 1
    else:
        return 0

def readLineRI(w, s):
    if s == 'MXE':
        line = {
        'ID': w[2],
        'GeneID': w[3],
        'geneSymbol': w[4],
        'chr': w[5],
        'strand': w[6],
        'SplicingStart': int(w[7]),
        'SplicingEnd': int(w[8]),
        'upstreamES': int(w[11]),
        'upstreamEE': int(w[12]),
        'downstreamES': int(w[13]),
        'downstreamEE': int(w[14]),
        'PValue': w[20],
        'FDR': w[21],
        'significant': w[29]}
    else:
        line = {
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

def read(filename, pG4, window, v, s, signi):
    RIpG4 = {}
    Densities = []
    cpt = 0
    with open(filename) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l and not l.startswith('newID2'):
                w = l.split('\t')
                w = readLineRI(w, s)
                if w['GeneID'] in pG4 and w['significant'] == signi:
                    cptG4 = 0
                    for G4 in pG4[ w['GeneID'] ]:
                        loc = getLocationpG4(w, pG4[ w['GeneID'] ][G4], window)
                        if loc != 0 :
                            cptG4 += 1
                            RIpG4[cpt] = {}
                            RIpG4[cpt].update(w)
                            RIpG4[cpt].update( pG4[ w['GeneID'] ][G4])
                            RIpG4[cpt]['Location'] = loc
                            cpt += 1
                    Densities.append(float(cptG4) / float( w['downstreamES'] +100 - w['upstreamEE'] -100) *1000)
    return RIpG4, Densities

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

def getpG4NearRI(path, window):
    pG4All = path + 'HS_pG4.csv'
    files = {'SE' : {'kunvRI' : path + 'SE/KUNV_SE.csv',
                    'sinvRI' : path + 'SE/SINV_SE.csv',
                    'zikvRI' : path + 'SE/ZIKV_SE.csv',
                    'yvfRI' : path + 'SE/YFV_SE.csv'},
            'MXE' : {'kunvRI' : path + 'MXE/KUNV_MXE.csv',
                            'sinvRI' : path + 'MXE/SINV_MXE.csv',
                            'zikvRI' : path + 'MXE/ZIKV_MXE.csv',
                            'yvfRI' : path + 'MXE/YFV_MXE.csv'},
            'A5SS' : {'kunvRI' : path + 'A5SS/KUNV_A5SS.csv',
                            'sinvRI' : path + 'A5SS/SINV_A5SS.csv',
                            'zikvRI' : path + 'A5SS/ZIKV_A5SS.csv',
                            'yvfRI' : path + 'A5SS/YFV_A5SS.csv'},
            'A3SS' : {'kunvRI' : path + 'A3SS/KUNV_A3SS.csv',
                            'sinvRI' : path + 'A3SS/SINV_A3SS.csv',
                            'zikvRI' : path + 'A3SS/ZIKV_A3SS.csv',
                            'yvfRI' : path + 'A3SS/YFV_A3SS.csv'},
            'RI' : {'kunvRI' : path + 'RI/KUNV_RI.csv',
                            'sinvRI' : path + 'RI/SINV_RI.csv',
                            'zikvRI' : path + 'RI/ZIKV_RI.csv',
                            'yvfRI' : path + 'RI/YFV_RI.csv'}}
    pG4 = importpG4(pG4All)
    for s in files:
        print(s)
        for v in files[s]:
            statD = {}
            # 1 = significant, 0 = non_significant
            RIpG4, densities = read(files[s][v], pG4, window, v, s, '1')
            print('\t'+v)
            # print(densities)
            densities = np.array(densities)
            print('\t\tMean = ', str(densities.mean()))
            print('\t\tÉcart type = ', str(densities.std()))
            RIpG4, densities = read(files[s][v], pG4, window, v, s, '0')
            densities = np.array(densities)
            print('\t\tMean = ', str(densities.mean()))
            print('\t\tÉcart type = ', str(densities.std()))



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
