#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
Copyright Université de Sherbrooke, Département de Biochimie et Département d'Informatique

contact: sarah.belhamiti@usherbrooke.ca




---------------------------------------------------

``ExtractInfoSubset`` **module description**:

From clustering information, retrieve group and sequences 
.. moduleauthor:: Sarah.Belhamiti

December 2017

"""

import csv, math, numpy
import string
import subprocess
import Bio.Align.Applications
import sys
import Bio
import os, shutil
import pylab
import scipy.sparse as sparse
import re
import random
import matplotlib as mpl 
import matplotlib.pyplot as plt
import numpy as np
import argparse
import numpy as np 
import math 
from Bio import SeqIO
from scipy import stats
from pylab import *
import __future__
from math import sqrt
from numpy import mean, std
import pandas as pd
import matplotlib.gridspec as gridspec
###########################################################################################################################################

def CreateDicoByRegionG4(filename, condition, dico):
	""" Dictionary of region G4, where for each region G4, we keep the information about the biotype, the localisation or the sequence
	    Parameters
	    ----------
	    filename : string
		name of file with all informations of G4 in the TRANSCRIPTOME
	    condition : string
		condition of information (the biotype, the localisation or the sequence)
	    dico : dictionary
		Dictionary of region G4

	    Returns
	    ----------
	    dico : dictionary
		Dictionary of region G4
	"""
	dico={}
	inputfile= open(filename,"r") # file opening for reading
	for line in inputfile: # for each transcript
		if ('InfoG4ByTranscript' not in line):
			words=line.split('\t')
			InfoG4ByTranscript=words[0].rstrip() 
			transcriptId=InfoG4ByTranscript.split('|')[0]
			chromosome=InfoG4ByTranscript.split('|')[1].split(':')[0]
			startG4=InfoG4ByTranscript.split('|')[1].split(':')[1].split('-')[0]
			endG4=InfoG4ByTranscript.split('|')[1].split(':')[1].split('-')[1]
			strand=InfoG4ByTranscript.split('|')[2]
			cGcC=words[1].rstrip() 
			G4Hunter=words[2].rstrip() 
			sequenceG4=words[3].rstrip() 
			G4NN=words[4].rstrip() 
			localisation=words[5].rstrip() 
			transcriptBiotype=words[6].rstrip() 
			InfoG4=chromosome+':'+startG4+'-'+endG4+'|'+strand
			if (dico.has_key(InfoG4) == False): # if dico not contain info for this region G4 and fot this parameters
				listParameter=[]
			else: # if dico already contain info for this region G4 and fot this parameters
				listParameter=dico.get(InfoG4)#get info already known
			if (condition=='biotype'):
				if (transcriptBiotype not in listParameter):
					listParameter.append(transcriptBiotype)
			elif (condition=='localisation'):
				if (localisation not in listParameter):
					listParameter.append(localisation)
			elif (condition=='sequenceG4'):
				if (sequenceG4 not in listParameter):
					listParameter.append(sequenceG4)
			dico[InfoG4]=listParameter
	return dico

######################################################################################################################################################
def CreateDicoByCluster(path, kM, dico):
	""" Dictionary of region G4 contain in each cluster
	    Parameters
	    ----------
	    path : string
		path of directory which contain file for ONE cluster contain in kM
	    kM : integer
		number max of cluster tested
	    dico : dictionary
		Dictionary of cluster

	    Returns
	    ----------
	    dico : dictionary
		Dictionary of cluster
	"""
	dico={}
	for filename in os.listdir(path+'/cluster_'+str(kM) ): # for each file of the directory
		numberCluster=filename.split("_")[1]
		if ('LIST' in filename):
			inputfile= open(path+'/cluster_'+str(kM)+'/'+filename,"r") # file opening for reading
			for line in inputfile: # for each transcript
				words=line.split(';')
				idRegion=words[0].rstrip() 
				if (dico.has_key(numberCluster) == False): # if dico not contain info for this region G4 and fot this parameters
					listRegion=[]
				else: # if dico already contain info for this region G4 and fot this parameters
					listRegion=dico.get(numberCluster)#get info already known
				if (idRegion not in listRegion):
						listRegion.append(idRegion)
				dico[numberCluster]=listRegion
	return dico
######################################################################################################################################################
def RetrievePosInList(liste, element):
	""" Give position of a element in a list
	    Parameters
	    ----------
	    list : liste
		list to test
	    element : string
		element searched in the list
	    Returns
	    ----------
	    position : integer
		position of the element in the list
	"""
	position=liste.index(element)
	return position

######################################################################################################################################################
def ExtractfilePercent(path, localisationName, RegionByCluster, LocalisationByRegion, kM):
	
	outputfile= open(path+'/cluster_'+str(kM)+'/cluster_TEST_correlation.csv',"w") # file opening for readingprint '\t'.join(str(e) for e in localisationName)
	header="cluster\t"+"\t".join(str(e) for e in localisationName)
	outputfile.write(header+"\n")
	sizeClusterAll=0
	for numberCluster, listIdRegion in RegionByCluster.items(): # for each cluster
		localisationCount=[]
		percent=[]
		for i in range(len(localisationName)):
			localisationCount.append(0)
			percent.append(0)
		sizeCluster=len(listIdRegion) # size of cluster = nbr of region G4 in it
		for idRegion in listIdRegion:
			listLocalisation=LocalisationByRegion.get(idRegion)
			for localisation in listLocalisation:
				pos=RetrievePosInList(localisationName, localisation) # position of localisation
				localisationCount[pos]=localisationCount[pos]+1# add +1 for this position and put it in the list
		for i in range(len(percent)):
			percent[i]=round(localisationCount[i]/float(sizeCluster)*100,2)
		sizeClusterAll=sizeClusterAll+sizeCluster
		line=str(numberCluster)+"\t"+"\t".join(str(e) for e in percent)
		outputfile.write(line+"\n")
	print "For k-means with k="+str(kM)+ " numbers regions = "+str(sizeClusterAll)
######################################################################################################################################################
def CalculMean(uneListe): 
	""" Give mean of one liste
	    Parameters
	    ----------
	    uneListe : list
		list of value
	    Returns
	    ----------
	     mean : float
		 mean of a list
	"""
	somme=0
	mean=0.0
	for element in uneListe:
		somme=somme+float(element)
	if (len(uneListe) !=0):
		mean=somme/float(len(uneListe))
	return mean
######################################################################################################################################################
def CalculVariance(uneListe, moy):
	""" Give variance of one liste, with a mean
	    Parameters
	    ----------
	    uneListe : list
		list of value
	     mean : float
		 mean of a list
	    Returns
	    ----------
	     var : float
		 variance of a liste
	"""
	variance=0.0
	for element in uneListe:
		variance=variance+((float(element)-moy)**2)
	if (len(uneListe) !=0):
		variance = variance/float(len(uneListe))
	return variance
######################################################################################################################################################
def CalculSquareroot(variance):
	""" Give stantard deviation from a variance
	    Parameters
	    ----------
	     var : float
		 variance of a liste
	    Returns
	    ----------
	     squareRoot : float
		 stantard deviation
	"""
	squareRoot=sqrt(variance) 
	return squareRoot
######################################################################################################################################################
def CalculStandardError(squareRoot, uneListe):	
	""" Give stantard error from a stantard deviation and a list which contain the values
	    Parameters
	    ----------
	     squareRoot : float
		 stantard deviation
	    uneListe : list
		list of value
	    Returns
	    ----------
	     standardError : float
		 stantard error
	"""
	if (uneListe):
		lenght=len(uneListe)
		standardError=squareRoot/sqrt(lenght) 
	else:
		standardError=0.0
	return standardError	

###########################################################################################################################################
def CreateFileSequence (path, kM,RegionByCluster, SequenceByRegion):
	""" for each cluster, create file of sequence fasta contained in this cluster
	    Parameters
	    ----------
	    path : string
		path of directory which contain file for ONE cluster contain in kM
	    kM : integer
		number max of cluster tested
	    RegionByCluster : dictionary
		Dictionary of cluster
	    Returns
	    ----------
	     standardError : float
		 stantard error
	"""
	for key, value in RegionByCluster.items(): # for each cluster
		inputfile= open(path+'/cluster_'+str(kM)+'/Sequence_C'+str(key)+".txt","w") # file opening for reading
		for idRegion in value:
			inputfile.write( ">"+ idRegion+"\n")
			inputfile.write( str("".join(SequenceByRegion.get(idRegion))) +"\n")

######################################################################################################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	parser.add_argument ('-p', '--path', default = '/home/local/USHERBROOKE/bels2814/Documents/These/clusterLocalisation/5UTR/Subset')
	parser.add_argument ('-kM', '--kMax', default = 2)
	parser.add_argument ('-f', '--file', default = '/home/local/USHERBROOKE/bels2814/Documents/These/HUMAN/all/HS_All_G4InTranscript.txt')
	
	return parser
	
######################################################################################################################################################
def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path	# directory which contain all the directory CHROMOSOME
	kM=arg.kMax
	filename=arg.file

	localisationName=['5','3','CDS','Intron','junction_5_CDS','junction_CDS_3','junction_CDS_CDS','junction_CDS_Intron','junction_Intron_CDS','ExonNC','IntronNC','junction_ExonNC_ExonNC','junction_ExonNC_IntronNC','junction_IntronNC_ExonNC']
	PercentByCluster={}
	

	RegionByCluster={}
	LocalisationByRegion={}
	SequenceByRegion={}
	

	RegionByCluster=CreateDicoByCluster(path, kM, RegionByCluster)	
	LocalisationByRegion=CreateDicoByRegionG4(filename, "localisation", LocalisationByRegion)
	SequenceByRegion=CreateDicoByRegionG4(filename, "sequenceG4", SequenceByRegion)
	CreateFileSequence (path,kM, RegionByCluster, SequenceByRegion)
	

	


	

	

main()
	
	

