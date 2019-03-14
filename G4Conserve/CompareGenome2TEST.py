#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
Copyright Université de Sherbrooke, Département de Biochimie et Département d'Informatique

contact: sarah.belhamiti@usherbrooke.ca

This software is a computer program whose compare region G4 found in transcriptome by G4screener to the OQs founds in genome.


---------------------------------------------------

``ConpareGenome`` **module description**:

The program create differents files output:

    * 
    * 

.. moduleauthor:: Sarah.Belhamiti

December 2017

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
	url = 'http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment='+chromosome+':'+start+','+end # connexion database
	response = urllib2.urlopen(url) # recuperation of informations
	fasta = response.read() # reading of information 
	return fasta


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
	segment=fasta.split('\n') # split all informations in the fonction 
	sequence=''
	for element in segment: # for each informations obtained
		if not (re.search('^<', element) or re.search('\n', element)): # if the information is the sequence (start with '>')
			sequence=sequence+element # get the information 
	if ('minus' in inputfile): # if on reverse strand, dna reverse --> sequence obtained transformed in dna complement
		dna = Seq(sequence, generic_dna)
		sequence_reverse= dna.reverse_complement()
		sequence=sequence_reverse
	return sequence


def ExtractInfo(inputfilename,dico):
	""" This fonction extract information from a file in a dictionnary, where key = chr:start-end|strand and value is sequence 
	    Parameters
	    ----------
	    inputfilename: string
		fonction which contain fasta sequence for each G4 or OQ
	    dico : dictionary
		
	    Returns
	    -------
	    dico : dictionary
	"""
	for line in inputfilename:	# for each OQs described in the file
		words=line.split('\t') # parsing by the separator, here '\t' (file bed)
		header=words[0].rstrip().replace("chr", "").replace("+1", "1")
		seq=words[1].rstrip()
		dico[header]=seq
	return dico

def CreationIndex(Index,headerTest,headerIndex):
	""" This fonction create an index for each chromosome and strand (key) , list of position of gene contained in this chromosome 
	    Parameters
	    ----------
	    Index : dictionary
		key= chr|strand, value = list of gene contained on this chromosome 
	    headerTest: string
		header of a gene, in format chr:start-end|strand
	    headerIndex: string
		header of a chromosome, in format chr|strand
	    Returns
	    -------
	    Index : dictionary
		key= chr|strand, value = list of gene contained on this chromosome 
	"""
	if (Index.has_key(headerIndex) == False): # if dictionary haven't this key 
		liste=[] # create empty list 
	else: # if dictionary have this key 
		liste=Index.get(headerIndex) # get the list already created
	coord=headerTest.split(':')[1].split('|')[0] # create coord of the gene contained in this chr
	liste.append(coord) # add this coord on the list
	Index[headerIndex]=liste
	return Index

def conditionInclusifGene (start, end, startBorder, endBorder):
	""" This fonction return a boolean about if a OQ is inclusif in a gene 
	    Parameters
	    ----------
	    start : string
		start OQ 
	    end : string
		end OQ 
	    startBorder : string
		start gene
	    endBorder : string
		end gene 
	    Returns
	    -------
	    condition : boolean
		reponse if inclusion or not 
	"""
	condition=False
	#geneExt=5000
	geneExt=0
	if (int(start) >= int(startBorder)-geneExt and int(end) <= int(endBorder)+geneExt):
		condition=True
	return condition


def conditionOverlap (startG4, endG4, startOQ, endOQ):
	""" This fonction return a boolean about if a OQ is overlapping with a G4 
	    Parameters
	    ----------
	    startG4 : string
		start G4 
	    endG4 : string
		end G4 
	    startOQ : string
		start OQ
	    endOQ : string
		end OQ 
	    Returns
	    -------
	    condition : boolean
		reponse if overlapping or not 
	"""
	test=False
	Incluse = (int(startOQ) >= int(startG4) and int(endOQ) <= int(endG4)) # if OQ inclus in G4r
	OverFilling = (int(startOQ) <= int(startG4) and int(endOQ) >= int(endG4)) # if OQ overfiling in G4r
	OverlapAmont = (int(startOQ) <= int(startG4) and int(endOQ) >= int(startG4) and int(endOQ) <= int(endG4)) # if OQ overlap in amont 
	OverlapAval = (int(startOQ) >= int(startG4) and int(startOQ) <= int(endG4) and int(endOQ) >= int(endG4)) # if OQ overlap in aval 
	
	if (Incluse or OverFilling or OverlapAmont or OverlapAval):
		test=True
	return test
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


	SeqOQ={}
	SeqG4={}
	IndexOQ={}
	IndexG4={}

	IndexGenome={}

	G4Ext=0
	#G4Ext=100
	geneExt=5000
	#geneExt=0

	condition='geneExt'
	print '\n'+ condition

	## creation liste of chromosome present for this specie
	listeChr=range(1,23)
	for i in range(len(listeChr)):
		listeChr[i]=str(listeChr[i])
	listeChr.append('X')
	listeChr.append('Y')



	## extraction information from G4 in genome and transcriptome = sequences
	inputfilename= open('/home/local/USHERBROOKE/bels2814/Documents/These/ChambersArticle/moreData/SeqOQ_2.txt',"r") # file opening for reading
	SeqOQ=ExtractInfo(inputfilename,SeqOQ)
	inputfilename= open('/home/local/USHERBROOKE/bels2814/Documents/These/ChambersArticle/moreData/SeqG4.txt',"r")	# file opening for reading
	SeqG4=ExtractInfo(inputfilename,SeqG4)

	## creation index for each G4 (OQ or G4r) 
	for headerOQ,seqOQ in SeqOQ.items():
		chrOQ=headerOQ.split(':')[0]
		strandOQ=headerOQ.split('|')[1]
		headerIndex=chrOQ+'|'+strandOQ
		if (chrOQ in listeChr):
			IndexOQ=CreationIndex(IndexOQ,headerOQ,headerIndex)
	for headerG4,seqG4 in SeqG4.items():
		chrG4=headerG4.split(':')[0]
		strandG4=headerG4.split('|')[1]
		headerIndex=chrG4+'|'+strandG4
		if (chrG4 in listeChr):
			IndexG4=CreationIndex(IndexG4,headerG4,headerIndex)


	## creation index for each chromosome
	inputfilename= open('/home/local/USHERBROOKE/bels2814/Documents/These/ChambersArticle/moreData/IndexGenome.txt',"r")	# file opening for reading
	for line in inputfilename:	# for each OQs described in the file
		words=line.split('\t') # parsing by the separator, here '\t' 
		idGene=words[0].rstrip()
		chrGene=words[1].rstrip()
		startGene=words[2].rstrip()
		endGene=words[3].rstrip()
		strandGene=words[4].rstrip()
		key=chrGene+'|'+strandGene
		if (chrGene in listeChr):
			if (IndexGenome.has_key(key) == False):
				liste=[]
			else:
				liste=IndexGenome.get(key)
			liste.append(startGene+'-'+endGene)
			IndexGenome[key]=liste


	

	## Get OQ in gene
	OQInGene=[]
	filename="/home/local/USHERBROOKE/bels2814/Documents/These/ChambersArticle/moreData/Test/OQInGene_"+condition+".txt"
	inputfile= open(filename,"r")	# file opening for reading
	for line in inputfile:	# for each line in the file (for each transcript)
		#chrm=line.rstrip().split(':')[0]
		#pos=line.rstrip().split(':')[1].split('|')[0]
		#strand=line.rstrip().split('|')[1]
		OQInGene.append(line.rstrip())


	## Get OQ in common
	OQCommon=[]
	filename="/home/local/USHERBROOKE/bels2814/Documents/These/ChambersArticle/moreData/Test/OQCommon_"+condition+".txt"
	inputfile= open(filename,"r")	# file opening for reading
	for line in inputfile:	# for each line in the file (for each transcript)
		OQCommon.append(line.rstrip())






	
	################### condition G4
	G4Common=[]
	for key, listeG4 in IndexG4.items():# for each chromosome and strand
		chrm=key.split("|")[0]
		strand=key.split("|")[1]
		for G4 in listeG4: 	# for each OQ
			commonG4=False
			startG4=int(G4.split('-')[0])-G4Ext
			endG4=int(G4.split('-')[1])+G4Ext
			listeOQ=sorted(IndexOQ.get(key)) # recup OQ in this chr and strand
			i=0
			while( i <len(listeOQ)  and commonG4 == False):
				startOQ=listeOQ[i].split('-')[0]
				endOQ=listeOQ[i].split('-')[1]
				if ((chrm+':'+startOQ+'-'+endOQ+'|'+strand) in OQInGene):# if OQ in gene:
					commonG4= conditionOverlap(startG4, endG4,startOQ, endOQ) # condition of OQ overlapp G4r
					if (commonG4 == True):
						if (G4 not in G4Common):
							G4Common.append(chrm+':'+G4+'|'+strand)
				i+=1


	outfilename= open('/home/local/USHERBROOKE/bels2814/Documents/These/ChambersArticle/moreData/Test/G4Common_'+condition+'.txt',"w")	# file opening for reading
	for element in G4Common:
		outfilename.write(element+'\n')	
	print 'extract G4 common done '
	print '----------------------------------------------------'

	print '\n'+ condition
	print '\nNumber OQ initial:	'+str(len(SeqOQ))
	print '\nNumber G4 initial:	'+str(len(SeqG4))
	print '\nNumber OQ genique:	'+str(len(OQInGene))
	print '\nNumber OQ overlap PG4r:	'+str(len(OQCommon))
	print '\nNumber PG4r overlap OQ:	'+str(len(G4Common))

			
main()
	
