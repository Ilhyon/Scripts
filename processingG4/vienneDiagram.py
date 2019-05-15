#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import argparse
import pandas as pd
import recurrentFunction as rF
import processToG4InTranscripts as merge

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'G4Annotation')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-sp', '--specie', default = \
        'yersinia_pestis_biovar_microtus_str_91001')
    parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
    parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
    parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
    parser.add_argument ('-E', '--EXTENSION', default = 100)
    parser.add_argument ('-W', '--WINDOW', default = 60)
    parser.add_argument ('-S', '--STEP', default = 10)
    return parser

def createDicoVenn():
    dico = {'G4NN' : {'g4NN' : 0.5,
                        'cGcC' : -4.5,
                        'g4H' : -4.5},
            'G4NN-G4H' : {'g4NN' : 0.5,
                            'cGcC' : -4.5,
                            'g4H' : 0.9},
            'G4NN-cGcC' : {'g4NN' : 0.5,
                            'cGcC' : 4.5,
                            'g4H' : -4.5},
            'G4H' : {'g4NN' : 0,
                    'cGcC' : -4.5,
                    'g4H' : 0.9},
            'G4H-cGcC' : {'g4NN' : 0,
                        'cGcC' : 4.5,
                        'g4H' : 0.9},
            'cGcC' : {'g4NN' : 0,
                    'cGcC' : 4.5,
                    'g4H' : -4.5},
            'cGcC-G4H-G4nn' : {'g4NN' : 0.5,
                                'cGcC' : 4.5,
                                'g4H' : 0.9}}
    return dico

def getStats(path, dicoParam):
    dicoTest = createDicoVenn()
    dicoRes = {}
    for subset in dicoTest:
        dicoParam.update(dicoTest[subset])
        dftmp = merge.mergeWindow(path, dicoParam)
        tmp = dftmp['Start'].map(str) +':'+ dftmp['End'].map(str)
        tmp = list(set(tmp))
        dicoRes.update({subset : tmp})
    return dicoRes

def getScores(path, dicoParam):
    dicoParam.update({'g4NN' : -5.0, 'cGcC' : -5.0, 'g4H' : -5.0})
    dicoRes = {}
    dftmp = merge.mergeWindow(path, dicoParam)
    dicoRes.update({'g4NN' : list(set(dftmp['G4NN'])),
                    'cGcC' : list(set(dftmp['cGcC'])),
                    'g4H' : list(set(dftmp['G4H']))})
    return dicoRes

def main(path, dicoParam, output):
    dico = getScores(path, dicoParam)
    df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in dico.items() ]))
    df.to_csv(path_or_buf=output, header=True, index=None, sep='\t', mode='a')

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    sp = arg.specie
    directory = path +'/'+ sp
    output = path  +'/All_sp/'+ sp + 'Score.csv'
    dicoParam = rF.createDicoParam(arg)
    main(directory, dicoParam, output)
