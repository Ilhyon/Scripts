#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import re
import os
import math
import random
import argparse
import pandas as pd
from Bio import SeqIO
from pprint import pprint
import Parser_gtf as pGTF
import recurrentFunction as rF

def writeFasta(fasta, outputDir, opt):
    """From a dictionary {id : seq}, write a fasta file.

    """
    if opt == 'Both':
        for type in fasta:
            output = open(outputDir + 'Sequences_' + type + '.fas', "w")
            for id in fasta[type]:
                output.write(id + "\n")
                nbLine = math.ceil( float( len(fasta[type][id]) ) / 60 )
                cpt1 = 0
                cpt2 = 60
                for i in range(0,int(nbLine)) :
                    output.write(fasta[type][id][cpt1:cpt2] + "\n")
                    # to have a new line after 60 characters
                    cpt1 += 60
                    cpt2 += 60
            output.close()
    else:
        output = open(outputDir + 'Sequences_' + opt + '.fas', "w")
        for id in fasta[type]:
            output.write(id + "\n")
            nbLine = math.ceil( float( len(fasta[type][id]) ) / 60 )
            cpt1 = 0
            cpt2 = 60
            for i in range(0,int(nbLine)) :
                output.write(fasta[type][id][cpt1:cpt2] + "\n")
                # to have a new line after 60 characters
                cpt1 += 60
                cpt2 += 60
        output.close()

def shuffleSeq(seq):
    """Shuffle a sequence.

    This function aims to shuffle a fasta sequence to randomize its sequence.
    The fasta sequence is imported from a fasta file, then converted into a
    list (one element corrresponds to one nucleotide), the list is shuffled and
    tehn joined with nothing to recreate the sequence.

    :param seq: sequence to shuffle.
    :type seq: string

    :returns: seq, sequence shuffled
    :rtype: string
    """
    seq = list(seq)
    random.shuffle(seq)
    seq = ''.join(seq)
    return seq

def reverseSequence(Sequence):
    """ Reverse complement a DNA sequence.

    :param Sequence: DNA sequence that will be reversed.
    :type Sequence: string

    :returns: Sequence, the initial DNA sequence but reverse complemented.
    :rtype: string
    """
    reverse = ""
    nucleotides = {'A' : 'T',
                    'T' : 'A',
                    'C' : 'G',
                    'G' : 'C'}
    for n in Sequence:
        if n in nucleotides:
            tmp = nucleotides[n]
        else :
            tmp = n # in some sequences there is many N or other letter
        reverse += tmp
    reverse = reverse[::-1]
    return reverse

def importFastaChromosome(filename, chr):
    """Import the fasta file containing an entire chromosome sequence.

    :param directory: directory name containing all chromosome fasta.
    :type directory: string
    :param chr: chromosome wanted.
    :type chr: string

    :returns: sequence, entire chromosome sequence.
    :rtype: string
    """
    with open(filename) as f: # file opening
        content = f.read()
        l = content.split('\n')
        if l[0].startswith('>'):
            header = l[0]
            sequence = "".join(l[1:])
    return sequence

def getDicoLocTrBt(dfChr):
    dicoLocTrBt = {}
    listTest = []
    listTest = dfChr['Chromosome'].astype(str) + '~' + \
        dfChr['Start'].astype(str) + '~' + \
        dfChr['End'].astype(str) + '~' + \
        dfChr['Strand'].astype(str) + '~' + \
        dfChr['Transcript'].astype(str) + '~' + \
        dfChr['Biotype'].astype(str) + '~' + \
        dfChr['Feature'].astype(str)
    for loc in listTest:
        w = loc.split('~')
        locId = w[0] + ':' + w[1] + '~' + w[2] + ':' + w[3]
        if locId not in dicoLocTrBt:
            dicoLocTrBt[locId] = []
        # add Tr + biotype
        dicoLocTrBt[locId].append(w[4] + '-' + w[5])
    return dicoLocTrBt

def createFasta(dfGTF, pathFasta, outputDir, opt):
    dicoFasta = { 'WT' : {}, 'Shuffled' : {} }
    fastaFiles = [f for f in os.listdir(pathFasta) if \
        os.path.isfile( os.path.join(pathFasta, f) )]
    chrList = list( set(dfGTF.Chromosome) )
    for chr in chrList:
        dfChr = pd.DataFrame()
        for file in fastaFiles:
            if ('chromosome.' + chr + '.fa' in file) or \
                ('chr' + chr + '.fa' in file) or \
                ('_group.' + chr + '.fa' in file):
                fastaFile = pathFasta + file
        chrSeq = importFastaChromosome(fastaFile, chr)
        dfChr = dfChr.append(dfGTF[ dfGTF.Chromosome == chr])
        dfChr['UniqID'] =  dfChr['Gene'].astype(str) + '~' + \
            dfChr['Start'].astype(str) + '~' + \
            dfChr['End'].astype(str) + '~' + \
            dfChr['Strand'].astype(str) + '~' + \
            dfChr['Feature'].astype(str)
        dicoLocTrBt = getDicoLocTrBt(dfChr)
        del dfChr['Transcript']
        del dfChr['Biotype']
        del dfChr['index1']
        del dfChr['Coords']
        del dfChr['Gene']
        dfChr = dfChr.reset_index(drop=True)
        dfChr = dfChr.drop_duplicates(subset=None, keep='first', inplace=False)
        dfChr = dfChr.reset_index(drop=True)
        dicoLoc = dfChr.set_index('UniqID').to_dict('index')
        for loc in dicoLoc:
            w = loc.split('~')
            locId = chr  + ':' + w[1] + '~' + w[2] + ':' + w[3]
            if dicoLoc[loc]['End'] != dicoLoc[loc]['Start']:
                dicoLoc[loc]['Start'] = int(dicoLoc[loc]['Start'])
                dicoLoc[loc]['End'] = int(dicoLoc[loc]['End'])
                if dicoLoc[loc]['Start'] > 0:
                    seq = chrSeq[ dicoLoc[loc]['Start'] - 1 : \
                        dicoLoc[loc]['End'] - 1 ]
                else:
                    negSeq = chrSeq[- -dicoLoc[loc]['Start'] - 1 :]
                    # sequence before the origin of replication
                    posSeq = chrSeq[ 0 : dicoLoc[loc]['End'] ]
                    # sequence after the origin of replication
                    seq = negSeq + posSeq
                if w[3] == 'junction':
                    seqUpstream = chrSeq[ dicoLoc[loc]['Start'] - 40 : \
                        dicoLoc[loc]['Start'] ]
                    seqDownstream = chrSeq[ dicoLoc[loc]['End'] : \
                        dicoLoc[loc]['End'] + 40]
                    seq = seqUpstream + seqDownstream
                if dicoLoc[loc]['Strand'] == '-1':
                    seq = reverseSequence(seq)
                randomSeq = shuffleSeq(seq)
                listTr = dicoLocTrBt[locId]
                id = '>' + w[0] +':'+ w[4] +':'+ locId +':'+ '|'.join(listTr)
                dicoFasta['WT'][id] = seq
                dicoFasta['Shuffled'][id] = seq
    writeFasta(dicoFasta, outputDir, opt)

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-o', '--option', default = 'Both')
    parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    sp = arg.specie
    path = arg.path
    opt = arg.option
    fastaFile = path + sp + '/Fasta/'
    outputDir= path + sp + '/'
    dfGTF = pGTF.importGTFdf(path + sp + '/' + sp + '.gtf', sp)
    createFasta(dfGTF, fastaFile, outputDir, opt)
