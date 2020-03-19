#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os
import argparse
from pprint import pprint

def writeOutputs(dicopG4, dicopG4Uniq, path):
    output = open(path + 'pG4RetainedIntronGeneName.csv', "w")
    for G4 in dicopG4:
        output.write(dicopG4[G4] + '\n')
    output.close()
    output = open(path + 'pG4RetainedIntronGeneName2.csv', "w")
    for G4 in dicopG4:
        l= dicopG4[G4].split('|')
        output.write( l[0] +','+ ('|').join(l[1:] )+ '\n')
    output.close()
    output = open(path + 'UniqpG4RetainedIntronGeneName.csv', "w")
    for G4 in dicopG4Uniq:
        output.write( (' ').join(dicopG4Uniq[G4]['gene']) + ',' + \
            dicopG4Uniq[G4]['line'] + '\n')
    output.close()

def readLinepG4(w):
    dicoL = {'idG4' : w[0],
            'cGcC' : w[1],
            'G4H' : w[2],
            'seq' : w[3],
            'G4NN' : w[4],
            'Location' : w[5],
            'Bt' : w[6]}
    return dicoL

def parsepG4Retained(file, dicoCorr):
    dicopG4 = {}
    dicopG4Uniq = {}
    with open(file) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines :
            if 'G4 uniq' not in l and l:
                w = l.split(',')
                l = readLinepG4(w)
                idTr = l['idG4'].split('|')[0]
                idG4 = ('|').join(l['idG4'].split('|')[1:])
                if idTr in dicoCorr:
                    if l['idG4'] not in dicopG4:
                        dicopG4[ l['idG4'] ] = dicoCorr[idTr] + ',' +  \
                                idG4 + ',' + \
                                l['cGcC'] + ',' + \
                                l['G4H'] + ',' + \
                                l['seq'] + ',' + \
                                l['G4NN'] + ',' + \
                                l['Location'] + ',' + \
                                l['Bt']
                    if idG4 not in dicopG4Uniq:
                        dicopG4Uniq[ idG4 ] =  {'gene' : [],
                                'line' : l['cGcC'] + ',' + l['G4H'] + ',' + \
                                l['seq'] + ',' + l['G4NN'] + ',' + \
                                l['Location'] + ',' + l['Bt']}
                    dicopG4Uniq[ idG4 ]['gene'].append(dicoCorr[idTr])
    return dicopG4, dicopG4Uniq

def readCorr(file):
    dicoCorr = {}
    with open(file) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if 'Gene stable ID' not in l and l:
                w = l.split('\t')
                if w[1] not in dicoCorr:
                    dicoCorr[ w[1] ] = w[2]
                else:
                    print(w[1])
    return dicoCorr

def main(path):
    correspondanceFile = path + 'GeneTrName.txt'
    pG4RetainedFile = path + 'pG4RetainedIntron.csv'
    dicoCorr = readCorr(correspondanceFile)
    dicopG4, dicopG4Uniq = parsepG4Retained(pG4RetainedFile, dicoCorr)
    writeOutputs(dicopG4, dicopG4Uniq, path)

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'getDataFig')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	main(path)
