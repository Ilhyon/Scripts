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

def main(path):
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
        print(v)
        print(len(dfpG4v[v]))




def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    main(path)
