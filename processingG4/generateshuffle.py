#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
    Copyright Universite of Sherbrooke, departement of biochemistry and
    departement    of computation.

Date:
    February 2020

Description:
    This script aims to generate all shufled fasta needed. It will generate one
    fasta file for each segment location but also for junction. The id of each
    fasta is like : biotype:Trid:Start-End:Strand, exept for intron where the
    biotype is at the end. To generate random sequences we retrieve the sequence
    of the location we want, we transformed it into a list, then we shufle it
    (change all position of all element randomly), to finally join the shuffled
    list with no characters to recreate the sequence.

Data availability:
    * fasta sequences used with G4RNA screener were download from Ensembl. The human assembly is GRCh38p12.
    * transcripts information were downloaded from ensembl via ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/

Abbreviation:
    * bt for biotype
    * tr for transcript
    * dico for dictionnary
    * w for word or words
    * dir for directory
    * cpt for compteur (counter)
"""

import os
import re
import math
import random
import argparse
from Bio import SeqIO
from pprint import pprint
import Parser_gtf as pGTF
import recurrentFunction as rF

def fromFasta(filename, outDir, ini, dico, dicoIntronGene):
    fastaOri = SeqIO.parse(open(filename), 'fasta')
    fastaShuf = {}
    for fasta in fastaOri:
        name, seq = fasta.id, str(fasta.seq)
        seq = shuffleSeq(seq)
        chr = name.split(':')[2]
        idLoc = chr+':'+name.split(':')[3] +'~'+name.split(':')[4]+ \
            ':'+name.split(':')[5]
        gene = dicoIntronGene[idLoc]
        for loc in dico[chr][gene]['intron']:
            if idLoc in loc:
                listTr = dico[chr][gene]['intron'][idLoc]
        name = '>'+gene+':junction:'+idLoc+':'+'|'.join(listTr)
        fastaShuf[name] = seq
    writeFasta(outputDir, fastaShuf, ini, 'junction')

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

def writeFasta(outDir, dicoFasta, ini, type):
    """From a dictionary {id : seq}, write a fasta file.

    :param outDir: name of the directory where the output file need to be writen.
    :type outDir: string
    :param dicoFasta: {id : seq}
    :type dicoFasta: dictionary
    :param chr: name of the chromosome.
    :type chr: string
    """
    output = open(outDir + 'Shuffle_' + ini + '_' + type + '.fas', "w")
    for id in dicoFasta:
        output.write(id + "\n")
        nbLine = math.ceil( float( len(dicoFasta[id]) ) / 60 )
        cpt1 = 0
        cpt2 = 60
        for i in range(0,int(nbLine)) :
            output.write(dicoFasta[id][cpt1:cpt2] + "\n")
            # to have a new line after 60 characters
            cpt1 += 60
            cpt2 += 60
    output.close()

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

def createFasta(dicoGTF, dicoIntron, pathFasta, ini, outputDir):
    """Create fasta file from the index of location.

    Each sequence corresponds to one location. The id is like :
    Gene:location:chr:start-end:Tr1-bt|Tr2-bt
    Location of only one nucleotides are not added to the fasta file.
    The sequence is retrieved from a fasta file containing the entire
    chromosome. This file have been dowloaded from the ensembl FTP.

    :param dico: aims to contain all location of a chromosome.
	:type dico: dictionary
    :param fastaFile: directory containing fasta of all chromosomes.
	:type fastaFile: string
    :param chr: name of the chromosome.
	:type chr: string
    :param outputDir: name of the output directory.
	:type outputDir: string
    """
    fastaRandom = {}
    fastaFiles = [f for f in os.listdir(pathFasta) if \
        os.path.isfile( os.path.join(pathFasta, f) )]
    for chr in dicoGTF:
        for file in fastaFiles:
            if ('chromosome.' + chr + '.fa' in file) or \
                ('chr' + chr + '.fa' in file) or \
                ('_group.' + chr + '.fa' in file):
                fastaFile = pathFasta + file
        chrSeq = importFastaChromosome(fastaFile, chr)
        for gene in dicoGTF[chr]:
            for locName in dicoGTF[chr][gene]:
                for location in dicoGTF[chr][gene][locName]:
                    # print(location)
                    start = int(location.split(':')[1].split('~')[0])
                    end = int(location.split(':')[1].split('~')[1])
                    if end != start:
                        strand = location.split(':')[2]
                        if start > 0:
                            seq = chrSeq[start-1:end-1]
                        else:
                            negSeq = chrSeq[- -start:]
                            # sequence before the origin of replication
                            posSeq = chrSeq[0:end]
                            # sequence after the origin of replication
                            seq = negSeq + posSeq
                        if strand == '-1':
                            seq = reverseSequence(seq)
                        randomSeq = shuffleSeq(seq)
                        listTr = dicoGTF[chr][gene][locName][location]
                        id = '>'+gene+':'+locName+':'+location+':'+'|'.join(listTr)
                        fastaRandom[id] = randomSeq
    writeFasta(outputDir, fastaRandom, ini, 'location')

def importIntron(filename):
    dico = {}
    with open(filename) as f:
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            words = l.split('\t')
            if len(words) > 1:
                idIntron = words[1]+':'+words[2]+'~'+words[3]+':'+words[4]
                if words[0] not in dico:
                    dico[ words[0] ] = [ ]
                dico[ words[0] ].append(idIntron)
    return dico

def addIntronToGTF(dicoGTF, dicoIntron, dicoTrGene):
    dicoIntronGene = {}
    for tr in dicoIntron:
        for intron in dicoIntron[tr]:
            gene = dicoTrGene[tr]['Gene']
            chr = intron.split(':')[0]
            strand = intron.split(':')[2]
            start = intron.split(':')[1].split('~')[0]
            end = intron.split(':')[1].split('~')[1]
            idLoc = chr+':'+start+'~'+end+':'+strand
            if 'intron' not in dicoGTF[chr][gene]:
                dicoGTF[chr][gene]['intron'] = { idLoc : []}
            if idLoc not in dicoGTF[chr][gene]['intron']:
                dicoGTF[chr][gene]['intron'][idLoc] = []
            dicoGTF[chr][gene]['intron'][idLoc].append(tr+'~'+dicoTrGene[tr]['Biotype'])
            dicoIntronGene.update({idLoc : gene})
    return dicoGTF, dicoIntronGene

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    sp = arg.specie
    path = arg.path
    ini = rF.setUpperLetter(sp)
    fastaFile = path+sp+'/Fasta/'
    fastaJunction = path+sp+'/'+ini+'_transcript_unspliced.txt'
    outputDir = path+sp+'/'
    dicoGTF = pGTF.importLocationFronGTF(path+sp+'/'+sp+'.gtf')
    dicoTrGene = pGTF.getDicoTrGene(path+sp+'/'+sp+'.gtf')
    dicoIntron = importIntron(path + sp +'/'+ ini + '_intron.txt')
    dicoLocation, dicoIntronGene = addIntronToGTF(dicoGTF, dicoIntron, dicoTrGene)
    createFasta(dicoGTF, dicoIntron, fastaFile, ini, outputDir)
    fromFasta(fastaJunction, outputDir, ini, dicoGTF, dicoIntronGene)
