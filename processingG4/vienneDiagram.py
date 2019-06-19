#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import getpG4
import argparse
import pandas as pd
import seaborn as sns
import recurrentFunction as rF
from matplotlib import pyplot as plt
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

def main(path, dicoParam, output, listSp):
    df = pd.DataFrame()
    for sp in listSp:
        shortSp = rF.getShortSpName(sp)
        ini = rF.setUpperLetter(sp)
        directory = path + sp
        dicoScore = merge.mergeWindow(directory, dicoParam, 'venn')
        for score in dicoScore:
            df = df.append(pd.DataFrame.from_dict({'ScoreValue' : dicoScore[score],
                                    'ScoreName' : [score] * len(dicoScore[score]),
                                    'Sp' : shortSp}))
    # ax = sns.boxplot(x='Sp', y='ScoreValue', data=df, palette=sns.color_palette("husl", 8))
    # ax = sns.swarmplot(x='Sp', y='ScoreValue', data=df, color="grey")
    # plt.show()
    fig, axs = plt.subplots(ncols=3)
    sns.boxplot(x='ScoreValue', y='Sp', data=df[df.ScoreName == 'G4NN'], ax=axs[0])
    sns.boxplot(x='ScoreValue', y='Sp', data=df[df.ScoreName == 'G4H'], ax=axs[1])
    sns.boxplot(x='ScoreValue', y='Sp', data=df[df.ScoreName == 'cGcC'], ax=axs[2])
    plt.xlim(0,40)
    # plt.xticks(rotation='vertical')
    plt.show()

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    sp = arg.specie
    listSp = ['yersinia_pestis_biovar_microtus_str_91001', 'saccharomyces_cerevisiae']
    directory = path +'/'+ sp
    output = path  +'/All_sp/'+ sp + 'Score.csv'
    dicoParam = rF.createDicoParam(arg)
    main(path, dicoParam, output, listSp)
