#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import getpG4
import argparse
import pandas as pd
import seaborn as sns
import recurrentFunction as rF
from matplotlib import pyplot as plt

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'Graphics')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
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

def main(path):
    df = pd.read_csv(path+'All_sp/Boxplot.csv', sep='\t')
    df.columns = ['ScoreValue', 'ScoreName', 'Sp']
    df = df[df.ScoreValue != 'ScoreValue'].dropna()
    df.ScoreValue = pd.to_numeric(df.ScoreValue)
    # data = df[df.ScoreName == 'G4H'].dropna()
    listSp = rF.createListSp()
    # sns.boxplot(x='ScoreValue', y='Sp', data=data, order=listSp)
    # plt.xlim(-0.1,2.5)
    # plt.draw()
    # plt.show()
    fig, axs = plt.subplots(ncols=3)
    sns.boxplot(x='ScoreValue', y='Sp', data=df[df.ScoreName == 'G4NN'],order=listSp, ax=axs[0])
    sns.boxplot(x='ScoreValue', y='Sp', data=df[df.ScoreName == 'G4H'],order=listSp, ax=axs[1])
    sns.boxplot(x='ScoreValue', y='Sp', data=df[df.ScoreName == 'cGcC'],order=listSp, ax=axs[2])
    # plt.xlim(0,40)
    plt.xscale('log')
    plt.draw()
    plt.show()
    fig.savefig(path+'All_sp/All_score.svg', dpi=500, format='svg', bbox_inches='tight')
    fig.savefig(path+'All_sp/All_score.png', dpi=500, format='png', bbox_inches='tight')

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    main(path)
