#!/usr/bin/env python
# -*- coding: utf-8 -*-:

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
    RIcoords = [ w['SplicingStart'] - window , w['SplicingEnd'] + window ]
    o = overlaps(pG4coords, RIcoords)
    if o > 0:
        if w['SplicingStart'] - window < w['upstreamES'] and w['SplicingStart'] + window > w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['downstreamEE'] ]
        elif w['SplicingStart'] - window < w['upstreamES'] and w['SplicingStart'] + window < w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['SplicingStart'] + window ]
        elif w['SplicingStart'] - window > w['upstreamES'] and w['SplicingStart'] + window > w['downstreamEE']:
            RIcoords = [ w['SplicingStart'] - window , w['downstreamEE'] ]
        else:
            RIcoords = [ w['SplicingStart'] - window , w['SplicingStart'] + window ]
        oUp = overlaps(pG4coords, RIcoords)
        if w['SplicingEnd'] - window < w['upstreamES'] and w['SplicingEnd'] + window > w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['downstreamEE'] ]
        elif w['SplicingEnd'] - window < w['upstreamES'] and w['SplicingEnd'] + window < w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['SplicingEnd'] + window ]
        elif w['SplicingEnd'] - window > w['upstreamES'] and w['SplicingEnd'] + window > w['downstreamEE']:
            RIcoords = [ w['SplicingEnd'] - window , w['downstreamEE'] ]
        else:
            RIcoords = [ w['SplicingEnd'] - window , w['SplicingEnd'] + window ]
        oDown = overlaps(pG4coords, RIcoords)
        if oUp > 0 and oDown <= 0:
            if w['strand'] == '+':
                return 'Up'
            else:
                return 'Down'
        elif oUp <= 0 and oDown > 0:
            if w['strand'] == '+':
                return 'Down'
            else:
                return 'Up'
        elif oUp > 0 and oDown > 0:
            return 'Both'
        else:
            return 0
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
    'SplicingStart': int(w[10]),
    'SplicingEnd': int(w[8]),
    'upstreamES': int(w[7]),
    'upstreamEE': int(w[10]),
    'downstreamES': int(w[11]),
    'downstreamEE': int(w[12]),
    'PValue': w[20],
    'FDR': w[21],
    'significant': w[29]}
    return line

def read(filename, pG4, window, v, signi):
    dicoRes = {'Up' : {'event' : [], 'pG4' : [], 'gene' : []},
        'Down' : {'event' : [], 'pG4' : [], 'gene' : []},
        'Both' : {'event' : [], 'pG4' : [], 'gene' : []}}
    output = []
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
                            loc = getLocationpG4(w, pG4[ w['GeneID'] ][G4], window)
                            if loc != 0:
                                # print(loc)
                                # print(v)
                                # print(w['geneSymbol'])
                                # print(w)
                                # print(pG4[ w['GeneID'] ][G4])
                                # print('----------')
                                dicoRes[loc]['event'].append(w['newID'])
                                dicoRes[loc]['gene'].append(w['geneSymbol'])
                                dicoRes[loc]['pG4'].append(str(pG4[ w['GeneID'] ][G4]['Start'])+'|'+\
                                    str(pG4[ w['GeneID'] ][G4]['End']))
                                output.append(loc+'\t'+w['geneSymbol']+'\t'+w['strand']+'\t'+\
                                    str(w['upstreamEE'])+'\t'+str(w['downstreamES'])+'\t'+\
                                    str(pG4[ w['GeneID'] ][G4]['Start'])+'\t'+str(pG4[ w['GeneID'] ][G4]['End']))

    return dicoRes, output

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
    files = {'kunv' : path + 'A5SS/KUNV_A5SS.csv',
            'sinv' : path + 'A5SS/SINV_A5SS.csv',
            'zikv' : path + 'A5SS/ZIKV_A5SS.csv',
            'yvf' : path + 'A5SS/YFV_A5SS.csv'}
    pG4 = importpG4(pG4All)
    for v in files:
        print(v)
        # 1 = significant, 0 = non_significant
        dicoRes, output = read(files[v], pG4, window, v, '0')
        # read(files[v], pG4, window, v, '0')
        outputF = open(path+v+'_A5SS0New.csv', "w")
        outputF.write( 'Location\tGeneSymbol\tStrand\tStartEvent\tEndEvent\tStartpG4\tEndpG4\n' )
        outputF.write( '\n'.join(output) )
        outputF.close()
        for loc in dicoRes:
            print('\t',loc)
            print('\t\tIl y a ', str(len(set(dicoRes[loc]['event']))), ' event avec au moins 1 pG4')
            print('\t\tIl y a ', str(len(set(dicoRes[loc]['gene']))), ' gene avec au moins 1 pG4')
            print('\t\tIl y a ', str(len(set(dicoRes[loc]['pG4']))), ' pG4')
        dicoRes, output = read(files[v], pG4, window, v, '1')
        # read(files[v], pG4, window, v, '0')
        outputF = open(path+v+'_A5SS1New.csv', "w")
        outputF.write( 'Location\tGeneSymbol\tStrand\tStartEvent\tEndEvent\tStartpG4\tEndpG4\n' )
        outputF.write( '\n'.join(output) )
        outputF.close()
        for loc in dicoRes:
            print('\t',loc)
            print('\t\tIl y a ', str(len(set(dicoRes[loc]['event']))), ' event avec au moins 1 pG4')
            print('\t\tIl y a ', str(len(set(dicoRes[loc]['gene']))), ' gene avec au moins 1 pG4')
            print('\t\tIl y a ', str(len(set(dicoRes[loc]['pG4']))), ' pG4')



def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-wDown', '--window', default = 100)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    window = arg.window
    getpG4NearRI(path, window)
