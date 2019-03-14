#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
Copyright Université de Sherbrooke, Département de Biochimie et Département d'Informatique

contact: sarah.belhamiti@usherbrooke.ca

This software is a computer program whose download sequence of OQ from the position.


---------------------------------------------------

``DownloadSeqOQ`` **module description**:

.. moduleauthor:: Sarah.Belhamiti


"""

import csv, math, numpy
import string
import subprocess
import Bio.Align.Applications
import sys
import Bio
import os
import scipy.sparse as sparse
import re
import random
import matplotlib as mpl 
import matplotlib.pyplot as plt
import numpy as np
import argparse
from Bio import SeqIO
from scipy import stats
from pylab import *
import numpy as np 
import urllib 
import urllib2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import gzip
from cogent.db.ensembl import HostAccount
from cogent.db.ensembl import Species
###########################################################################################################################################
def Fasta(chromosome,start,end):
	""" This fonction is to retrieve the fasta from chromosomique position for one sequence.Doesn't matter the sens of strand (always forward here )
	    Parameters
	    ----------
	    chromosome : string
		number of chromosome
	    start : string
		start sequence
	    end : string
		end sequence
	    Returns
	    -------
	    fasta: fonction
		fonction which contain fasta sequence from database
	"""
	url = 'http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment='+chromosome+':'+start+','+end
	response = urllib2.urlopen(url)
	fasta = response.read()
	return fasta
###########################################################################################################################################
def Sequence(fasta, inputfile):
	""" This fonction is to recreate the fasta from information doanloaded via database
	    Parameters
	    ----------
	    fasta: fonction
		fonction which contain fasta sequence from database
	    inputfile : string
		name of file which contain informations of positions
	    Returns
	    -------
	    sequence: string
		sequence obtained from the chromosomique position 
	"""
	segment=fasta.split('\n')
	sequence=''
	for element in segment:
		if not (re.search('^<', element) or re.search('\n', element)):
			sequence=sequence+element
	if ('minus' in inputfile): # if on reverce strand, dna reverse
		dna = Seq(sequence, generic_dna)
		sequence_reverse= dna.reverse_complement()
		sequence=sequence_reverse
	return sequence

######################################################################################################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	parser.add_argument ('-pOQ', '--pathOQ', default = '/home/local/USHERBROOKE/bels2814/Documents/These/ChambersArticle/moreData')

	return parser
######################################################################################################################################################
def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	pathOQ=arg.pathOQ	# directory which contain all the directory chromosome

	OqFasta={} # dictionnary of fasta squence for each Oq
	
	output=open(pathOQ+'/SeqOQ.txt',"w")
	fileList=['minus', 'plus'] # critearia of discrovered in forward or reverse strand
	for i in range(len(fileList)):
		criteria=fileList[i]
		inputfile=pathOQ+'/GSE63874_Na_K_'+criteria+'_hits_intersect.bed' # get name of file
	
		if (criteria == 'minus'):
			strand='-1'
		else:
			strand='+1'
		inputfilename= open(inputfile,"r")	# file opening for reading
		for line in inputfilename:	# for each OQs described in the file
			words=line.split('\t') # parsing by the separator, here '\t' (file bed)
			chromosomeOQ=words[0].rstrip()
			startOQ=words[1].rstrip()
			endOQ=words[2].rstrip()
			fasta=Fasta(chromosomeOQ,startOQ,endOQ) # get fasta file url
			sequence=Sequence(fasta, inputfile) # get sequence for this coordinates
			header=chromosomeOQ+':'+startOQ+'-'+endOQ+':'+strand
			OqFasta[header]=sequence
			output.write(header+'\t'+str(sequence).upper()+'\n')
main()
	
