#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import re
import os
import argparse
import ListrG4
from Bio import SeqIO
from pprint import pprint

def readFasta(fastaF):
    dicoG = {}
    G = "([G,g]*)"
    GG = "(G{2,}|g{2,})"
    canonical = r"(?i)(?P<run>g{3,})(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)"
    fasta = SeqIO.parse(open(fastaF),'fasta')
    for f in fasta:
        name, sequence = f.id, str(f.seq)
        listG = re.findall(GG, sequence)
        listG = list(filter(None, listG))
        if sum([len(x) for x in listG]) > 6 :
            dicoG[name] = listG
    pprint(dicoG)
    print('Number of rG4 with less then 4 GG :', len(dicoG))

def findNotMatching(fastaF):
    cpt = 0
    dicoG = {}
    matching = ListrG4.listrG4()
    fasta = SeqIO.parse(open(fastaF),'fasta')
    GG = r"(G{2,}|g{2,})"
    # minMotif = r"([G,g]{2,}[A-Z,a-z]{1,}){3,}"
    # minMotif = r"([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,}[A-Z,a-z]{1,})"
    minMotif = r"([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,})"
    for f in fasta:
        name, sequence = f.id, str(f.seq)
        if name not in matching:
            cpt += 1
            listG = re.findall(minMotif, sequence)
            if not listG:
                dicoG[name] = sequence
                print(sequence)
                print(listG)
                print('--------------')
            # listG = re.findall(minMotif, sequence)
            # print(sequence)
            # print(listG)
            # print('-----------------')
            # listG = list(filter(None, listG))
            # print(listG)
            # if sum([len(x) for x in listG]) <= 8:
            #     dicoG[name] = listG
            # else:
            #     print('>',name,'\n',sequence)
    # pprint(dicoG)
    print('Number of rG4 with less then 4 GG :', len(dicoG))
    print('Total not matching G4', cpt)

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'readBlast')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    fastaFile = path+'cdt_K.fas'
    findNotMatching(fastaFile)
