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
    RIcoords = [ w['Exon1S'] - window , w['Exon2E'] + window ]
    o = overlaps(pG4coords, RIcoords)
    if o > 0:
        if w['Exon1S'] - window < w['upstreamES'] and w['Exon1S'] + window > w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['downstreamEE'] ]
        elif w['Exon1S'] - window < w['upstreamES'] and w['Exon1S'] + window < w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['Exon1S'] + window ]
        elif w['Exon1S'] - window > w['upstreamES'] and w['Exon1S'] + window > w['downstreamEE']:
            RIcoords = [ w['Exon1S'] - window , w['downstreamEE'] ]
        else:
            RIcoords = [ w['Exon1S'] - window , w['Exon1S'] + window ]
        oUp1 = overlaps(pG4coords, RIcoords)
        if w['Exon1E'] - window < w['upstreamES'] and w['Exon1E'] + window > w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['downstreamEE'] ]
        elif w['Exon1E'] - window < w['upstreamES'] and w['Exon1E'] + window < w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['Exon1E'] + window ]
        elif w['Exon1E'] - window > w['upstreamES'] and w['Exon1E'] + window > w['downstreamEE']:
            RIcoords = [ w['Exon1E'] - window , w['downstreamEE'] ]
        else:
            RIcoords = [ w['Exon1E'] - window , w['Exon1E'] + window ]
        oDown1 = overlaps(pG4coords, RIcoords)
        if w['Exon2S'] - window < w['upstreamES'] and w['Exon2S'] + window > w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['downstreamEE'] ]
        elif w['Exon2S'] - window < w['upstreamES'] and w['Exon2S'] + window < w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['Exon2S'] + window ]
        elif w['Exon2S'] - window > w['upstreamES'] and w['Exon2S'] + window > w['downstreamEE']:
            RIcoords = [ w['Exon2S'] - window , w['downstreamEE'] ]
        else:
            RIcoords = [ w['Exon2S'] - window , w['Exon2S'] + window ]
        oUp2 = overlaps(pG4coords, RIcoords)
        if w['Exon2E'] - window < w['upstreamES'] and w['Exon2E'] + window > w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['downstreamEE'] ]
        elif w['Exon2E'] - window < w['upstreamES'] and w['Exon2E'] + window < w['downstreamEE']:
            RIcoords = [ w['upstreamES'] , w['Exon2E'] + window ]
        elif w['Exon2E'] - window > w['upstreamES'] and w['Exon2E'] + window > w['downstreamEE']:
            RIcoords = [ w['Exon2E'] - window , w['downstreamEE'] ]
        else:
            RIcoords = [ w['Exon2E'] - window , w['Exon2E'] + window ]
        oDown2 = overlaps(pG4coords, RIcoords)
        if oUp1 > 0 or oDown1 > 0 or oUp2 > 0 or oDown2 > 0:
            return 'hit'
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
    'Exon1S': int(w[7]),
    'Exon1E': int(w[8]),
    'Exon2S': int(w[9]),
    'Exon2E': int(w[10]),
    'upstreamES': int(w[11]),
    'upstreamEE': int(w[12]),
    'downstreamES': int(w[13]),
    'downstreamEE': int(w[14]),
    'PValue': w[20],
    'FDR': w[21],
    'significant': w[31]}
    return line

def read(filename, pG4, window, v, signi):
    dicoRes = {'event' : [], 'pG4' : [], 'gene' : []}
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
                                dicoRes['event'].append(w['newID'])
                                dicoRes['gene'].append(w['geneSymbol'])
                                dicoRes['pG4'].append(str(pG4[ w['GeneID'] ][G4]['Start'])+'|'+\
                                    str(pG4[ w['GeneID'] ][G4]['End']))
                                output.append(loc+'\t'+pG4[ w['GeneID'] ][G4]['loc']+'\t'+w['geneSymbol']+'\t'+w['strand']+'\t'+w['chr']+'\t'+\
                                    str(w['upstreamEE'])+'\t'+str(w['downstreamES'])+'\t'+\
                                    str(pG4[ w['GeneID'] ][G4]['Start'])+'\t'+str(pG4[ w['GeneID'] ][G4]['End'])+'\t'+\
                                    pG4[ w['GeneID'] ][G4]['Seq']+'\t'+\
                                    str(w['upstreamES'])+'\t'+str(w['upstreamEE'])+'\t'+\
                                    str(w['downstreamES'])+'\t'+str(w['downstreamEE'])+'\t'+\
                                    ' '.join(pG4[ w['GeneID'] ][G4]['Tr']))

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
                if w[9]+'|'+w[1] not in pG4[ w[9] ]:
                    pG4[ w[9] ][ w[9]+'|'+w[1] ] = {'Start' : int(start),
                        'End' : int(end),
                        'Chr' : w[1].split(':')[0],
                        'Strand' : w[2],
                        'cGcC' : w[3],
                        'G4H' : w[4],
                        'G4nn' : w[6],
                        'Seq' : w[5],
                        'Bt' : w[8],
                        'loc' : w[7],
                        'Tr' : [],}
                pG4[ w[9] ][ w[9]+'|'+w[1] ]['Tr'].append(w[0])
    return pG4

def getpG4NearRI(path, window):
    pG4All = path + 'HS_pG4.csv'
    files = {'kunvRI' : path + 'MXE/KUNV_MXE.csv',
            'sinvRI' : path + 'MXE/SINV_MXE.csv',
            'zikvRI' : path + 'MXE/ZIKV_MXE.csv',
            'yvfRI' : path + 'MXE/YFV_MXE.csv'}
    pG4 = importpG4(pG4All)
    for v in files:
        print(v)
        # 1 = significant, 0 = non_significant
        dicoRes, output = read(files[v], pG4, window, v, '0')
        # read(files[v], pG4, window, v, '0')
        outputF = open(path+v+'_MXE_NonSigni.csv', "w")
        outputF.write( 'Location\tpG4Location\tGeneSymbol\tStrand\tchr\tStartEvent\tEndEvent\tStartpG4\tEndpG4\tpG4Sequence\tE1S\tE1E\tE2S\tE2E\tTranscripts\n' )
        outputF.write( '\n'.join(output) )
        outputF.close()
        print('\t\tIl y a ', str(len(set(dicoRes['event']))), ' event avec au moins 1 pG4')
        print('\t\tIl y a ', str(len(set(dicoRes['gene']))), ' gene avec au moins 1 pG4')
        print('\t\tIl y a ', str(len(set(dicoRes['pG4']))), ' pG4')
        dicoRes, output = read(files[v], pG4, window, v, '1')
        # read(files[v], pG4, window, v, '0')
        outputF = open(path+v+'_MXE_Signi.csv', "w")
        outputF.write( 'Location\tpG4Location\tGeneSymbol\tStrand\tchr\tStartEvent\tEndEvent\tStartpG4\tEndpG4\tpG4Sequence\tE1S\tE1E\tE2S\tE2E\tTranscripts\n' )
        outputF.write( '\n'.join(output) )
        outputF.close()
        print('\t\tIl y a ', str(len(set(dicoRes['event']))), ' event avec au moins 1 pG4')
        print('\t\tIl y a ', str(len(set(dicoRes['gene']))), ' gene avec au moins 1 pG4')
        print('\t\tIl y a ', str(len(set(dicoRes['pG4']))), ' pG4')



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
