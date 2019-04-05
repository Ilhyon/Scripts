#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import Parser_gtf as pars
from Bio import SeqIO
from pprint import pprint

def getUTR(dicoTr, dicoProt):
    for tr in dicoTr:
        if tr in dicoProt:
            if dicoTr[tr]['Start'] == dicoProt[tr]['Start'] or \
                dicoTr[tr]['End'] == dicoProt[tr]['End']:
                lengthTr = abs(int(dicoTr[tr]['Start']) - \
                                int(dicoTr[tr]['End']))
                diffLength = lengthTr - len(dicoProt[tr]['Sequence']) +1
                if diffLength >= 1 :
                    truc = 0
                    # if '5UTR' in dicoTr[tr] or '3UTR' in dicoTr[tr]:
                    #     print 'Patate'
                elif diffLength == 0:
                    if '5UTR' in dicoTr[tr] or '3UTR' in dicoTr[tr]:
                        pprint(dicoTr[tr])
                        print dicoProt[tr]['Sequence']
                        print '-------------'
                # else:
                    # print 'whaaat'

def getInfoTr(description):
    description = description.split(' ')
    dicoElem = {}
    for elem in description:
        if re.search("chromosome:", elem) or \
            re.search("plasmid:", elem):
            info = elem.split(':')
            dicoElem = {'Chromosome' : info[2],
                        'Start' : info[3],
                        'End' : info[4],
                        'Strand' : info[5]}
    return dicoElem

def importProteinFasta(filename):
    dicoProt = {}
    for record in SeqIO.parse(filename, "fasta"):
        dicoProt[record.id] = {'Sequence' : record.seq}
        dicoProt[record.id].update(getInfoTr(record.description))
    return dicoProt

def main():
    filename = '/home/anais/Documents/Data/Genomes/' + \
                'drosophila_melanogaster/drosophila_melanogaster_cds.fa'
    dicoTr, dicoGene = pars.importGTF('/home/anais/Documents/Data/Genomes/' + \
                        'drosophila_melanogaster/drosophila_melanogaster.gtf')
    del dicoGene
    dicoProt = importProteinFasta(filename)
    getUTR(dicoTr, dicoProt)

if __name__ == '__main__':
	main()
