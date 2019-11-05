#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import sys
import os
import re
import argparse
from Bio.Blast import NCBIXML

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'readBlast')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-nHit', '--nHit', default = 5)
    return parser

def GetBestHit(path, nHit):
    filelist = sorted(os.listdir(path))# recuperation of all file name in the directory
    bestHit = {} # initialisation of the dictionary
    cptEmpty = 0
    for blastFile in filelist:
        tmp = blastFile.split(':')
        chrQ = tmp[2].split('r')[1]
        strandQ = tmp[5].rstrip()
        if strandQ == '+':
            strandQ = '1'
        else:
            strandQ = '-1'
        namequery = blastFile.split('.')[0]
        blastFile = path+blastFile
        blast = NCBIXML.parse(open(blastFile,'rU'))
        for record in blast:
            if record.alignments:
                hit = False
                hitN = 0
                if len(record.alignments) < nHit:
                    maxHit = len(record.alignments) -1
                else:
                    maxHit = nHit
                while not hit and hitN < maxHit:
                    # print(record.alignments[0])
                    strandT = record.alignments[hitN].title.split('|')[3]
                    chrT = record.alignments[hitN].title.split('|')[2].split(':')[0].split(' ')[1]
                    hitN += 1
                    if strandT == strandQ and chrT == chrQ:
                        hit = True
                        print('Ah')
                        # print(blastFile)
                        # print(record.alignments[hitN])
                # print(record.alignments[0].title)
                # print(record.alignments[0].hsps[0].expect)
            else:
                cptEmpty += 1
    print(cptEmpty)



if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    nHit = int(arg.nHit)
    GetBestHit(path, nHit)
