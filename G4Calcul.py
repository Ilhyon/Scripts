#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
Copyright Université de Sherbrooke, Département de Biochimie et Département d'Informatique

contact: sarah.belhamiti@usherbrooke.ca

This software is a computer program whose calcul the G4 region of one CHROMOSOME (or all) of one specie.


---------------------------------------------------

``G4Calcul`` **module description**:

From a tabulate files which contain G4 regions detected in one CHROMOSOME (or all chromosomes), this module will study the the distribution of the regions in the transcriptome. 
The program create differents files output:

    * 
    * 

.. moduleauthor:: Sarah.Belhamiti

December 2017

"""

import csv, math, numpy
import string
import subprocess
#import Bio.Align.Applications
import sys
#import Bio
import os, shutil
#import pylab
import scipy.sparse as sparse
import re
import random
import matplotlib as mpl 
import matplotlib.pyplot as plt
import numpy as np
import argparse
import numpy as np 
import math 
#from Bio import SeqIO
from scipy import stats
from pylab import *
import __future__
from math import sqrt
from numpy import mean, std
#import pandas as pd
#import matplotlib.gridspec as gridspec
######################################################################################################################################################
def GetLengthFraction(positionA,positionB):
	""" This fonction is to define the lenght between two positions (positionA and positionB).Doesn't matter the order of position (positionA can be > positionB)
	    Parameters
	    ----------
	    positionA : string
		position start for the position for exemple
	    positionB : string
		position end for the position for exemple
	    Returns
	    -------
	    length: integer
		lenght between two positions
	"""
	length=0
	if ((positionA and positionB) != ""): # because some transcript doesn't have 5' for exemple (and start5 and end5 will be == '')
		length=(abs(int(positionA)-int(positionB)))
	return int(length)
######################################################################################################################################################
def GetLengthCodant(exonList):
	""" This fonction is to define the lenght of exonique part.Can be use for the intronique part.
	    Parameters
	    ----------
	    exonList : list
		list of couple exon contain in one transcript
	    Returns
	    -------
	    lengthCodant: integer
		lenght of exonique part
	"""
	lengthCodant=0
	if ('' not in exonList): # if transcript contain exon
		for coupleExon in exonList: # for each exon
			startExon=coupleExon.split('-')[0] 
			endExon=coupleExon.split('-')[1]
			lengthCodant+=GetLengthFraction(startExon,endExon)# add the lenght of this exon to the total (lengthCodant)
	return lengthCodant 
######################################################################################################################################################
def GetLengthTotal(exonList):
	""" This fonction is to define the lenght total of a transcript based on the list of exon
	    Parameters
	    ----------
	    exonList : list
		list of couple exon contain in one transcript
	    Returns
	    -------
	    lengthCodant: integer
		lenght of exonique part
	"""
	lengthTotal=0
	allExonsBorder=[] # all borders of exons contained in a transcript
	for coupleExon in exonList:
		allExonsBorder.append(int(coupleExon.split('-')[0]))
		allExonsBorder.append(int(coupleExon.split('-')[1]))
	lengthTotal=GetLengthFraction(min(allExonsBorder),max(allExonsBorder))
	return lengthTotal
######################################################################################################################################################
def GetPositionForLocalisation(localisation, typeFamily):
	""" This fonction is to give the number of position for the localisation depending to the type of trasnscript (coding or not)
		See the domentation of GetLenghtByTranscript() to know the corresponding of those positions
	    Parameters
	    ----------
	    localisation : string
		localisation of the G4 on the transcript
	    typeFamily : string
		transcript biotype
	    Returns
	    -------
	    position: integer
		position of G4 depending to the localisation
	"""
	if (typeFamily =='Coding'):
		if(localisation == '5'):
			position=0
		elif(localisation == 'CDS'):
			position=1
		elif(localisation == '3'):
			position=2
		elif(localisation == 'Intron'):
			position=4
		elif(localisation == 'junction_5_CDS'):
			position=5
		elif(localisation == 'junction_CDS_Intron'): 
			position=6
		elif(localisation == 'junction_Intron_CDS'):
			position=7
		elif(localisation == 'junction_CDS_3'):
			position=8
		elif(localisation == 'junction_CDS_CDS'):
			position=10
		#### position 3 for total exon (5, CDS, 3)
		#### position 9 for total in transcript entire (all except junction_CDS_CDS)
	else: # number G4 in exon, intron, junction CDS, junction i, total, junction cds_cds
		if(localisation == 'ExonNC'):
			position=0
		elif(localisation == 'IntronNC'):
			position=1
		elif(localisation == 'junction_ExonNC_IntronNC'): 
			position=2
		elif(localisation == 'junction_IntronNC_ExonNC'):
			position=3
		elif(localisation == 'junction_ExonNC_ExonNC'):
			position=5
		#### position 4 for total in transcript entire (all except junction_CDS_CDS))
	#~ print(position)
	return position
		
	
######################################################################################################################################################
def GetAnnotationTranscript(filename,Coding,BiotypeByTranscript):
	""" This fonction is to define which transcript contain a good annotation or not from Ensembl
	    Creation of a dico, where key = id of the transcript and value= bolean value (False if the annotation isn't good)
	    Parameters
	    ----------
	    filename : string
		name of the file which contain informatiob fr each transcript (except the transcript's biotype)
	    Coding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    BiotypeByTranscript : dico
		dictionary, where key =id of transcript and value = biotype of this transcript
	    Returns
	    -------
	    AnnotationTranscript: dico
		dictionary with information of annotation for the transcripts By Ensembl (True if good anotation, else False)
	"""
	AnnotationTranscript={} # creation dico empty
	inputfile= open(filename,"r") ## file opening with informations of transcipts
	for line in inputfile: # for each transcript 
		words=line.split('|')
		transcriptId=words[0].rstrip() 
		start5=words[7].rstrip()
		end5=words[8].rstrip() 
		start3=words[9].rstrip() 
		end3=words[10].rstrip() 
		answer=True # variable answer by defaul True
		transcriptBiotype= ''
		transcriptBiotype=BiotypeByTranscript.get(transcriptId)	## get transcriptBiotype of this TRanscriptId
		if (transcriptBiotype not in Coding): # if transcript not a protein coding
			if (start5!='' or end5!='' or start3!='' or end3!='' ): # but if transcript has a 5'UTR or an 3' UTR
				answer=False # has not a good annotation in Ensembl	
		AnnotationTranscript[transcriptId]=answer
	inputfile.close()
	return AnnotationTranscript
######################################################################################################################################################
def CreateDictionaryBiotypeByTranscript (filename):
	""" Create dictionary with as value the biotype of the transcript for all transcripts of the chromosome

	    Parameters
	    ----------
	    filename : string
		name of file which contain the liste of biotype for each transcript of this chromosome

	    Returns
	    -------
	    dico
		dictionary
		where key = ID of transcript, value = biotype of this transcript
	"""
	dico={}	#Creation dictionary empty, where will be : key = transcript, value = biotype of transcript
	inputfile= open(filename,"r")	# file opening for reading
	for line in inputfile:	# for each line in the file (for each transcript)
		words=line.split('|')	# parsing by the separator, here '|'
		#~ print(line)
		transcript=words[1].rstrip()	#I D of transcript
		transcriptBiotype=words[3].rstrip() # biotype of the transcript
		#if (dico.has_key(transcript) == False):	# if transcript not contain his biotype in the dico
		dico[transcript]=transcriptBiotype # create an entry (biotype) for this transcript
	inputfile.close()
	return dico	# return dico with all transcript of this chromosome
######################################################################################################################################################
def GetLenghtByTranscript(filename, BiotypeByTranscript,AnnotationTranscript,Coding ):
	""" Create dictionary of lenght for each transcript
	if transcript is coding ,format =  [length 5,length cds,length 3,lenght exon,lenght intron,number junction 5_codant,number junction codant_Intron,number junction Intron_codant,number junction codant_3,length total, number junction CDS_CDS ]
	else format = [lenght exon,lenght intron,number junction codant_Intron,number junction Intron_codant,length total, number junction CDS_CDS ]
	    Parameters
	    ----------
	    filename : string
		name of the file which contain information for each transcript (except the transcript's biotype)
 	    BiotypeByTranscript : dico
		dictionary, where key =id of transcript and value = biotype of this transcript
	    AnnotationTranscript: dico
		dictionary with information of annotation for the transcripts By Ensembl (True if good anotation, else False)
	    Coding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    Returns
	    -------
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	"""
	LenghtByTranscript={}
	inputfile= open(filename,"r") # file opening for reading
	for line in inputfile: # for each transcript
			words=line.split('|')
			transcriptId=words[0].rstrip() 
			geneId=words[1].rstrip() 
			chromosome=words[2].rstrip() 
			strand=words[3].rstrip() 
			biotypeGene=words[4].rstrip() 
			exonList=words[5].rstrip().split(";") ## transform in array
			intronList=words[6].rstrip().split(";") ## transform in array
			start5=words[7].rstrip()
			end5=words[8].rstrip()
			start3=words[9].rstrip()
			end3=words[10].rstrip()
			annotationTranscript=AnnotationTranscript.get(transcriptId)
			transcriptBiotype=BiotypeByTranscript.get(transcriptId) # get biotype of this transcript
			if (annotationTranscript == True and transcriptBiotype is not None): # if annotation of this transcript is correct
				if (transcriptBiotype in Coding):
					infoTranscript=[0,0,0,0,0,0,0,0,0,0,0] #format [length 5,length cds,length 3,lenght exon,lenght intron,number junction 5_codant,number junction codant_Intron,number junction Intron_codant,number junction codant_3,length total, number junction CDS_CDS ]
					infoTranscript[0]=GetLengthFraction(start5,end5) #length 5
					infoTranscript[1]=GetLengthCodant(exonList)-GetLengthFraction(start5,end5)-GetLengthFraction(start3,end3) #length cds
					infoTranscript[2]=GetLengthFraction(start3,end3) #length 3
					infoTranscript[3]=GetLengthCodant(exonList) #lenght exon
					infoTranscript[4]=GetLengthCodant(intronList) #lenght intron
					infoTranscript[5]=1 # number junction 5_codant
					infoTranscript[6]=len(exonList)-1 #  number junction codant_Intron
					infoTranscript[7]=len(exonList)-1 #number junction Intron_codant
					infoTranscript[8]=1 # number junction codant_3
					infoTranscript[9]=GetLengthTotal(exonList) #length total
					infoTranscript[10]=len(exonList)-1#number junction CDS_CDS
					LenghtByTranscript[transcriptId+'-'+transcriptBiotype]=infoTranscript
				else:
					infoTranscript=[0,0,0,0,0,0] #format [lenght exon,lenght intron,number junction codant_Intron,number junction Intron_codant,length total, number junction CDS_CDS ]
					infoTranscript[0]=GetLengthCodant(exonList) #lenght exon
					infoTranscript[1]=GetLengthCodant(intronList) #lenght intron
					infoTranscript[2]=len(exonList)-1 #  number junction codant_Intron
					infoTranscript[3]=len(exonList)-1 #number junction Intron_codant
					infoTranscript[4]=GetLengthTotal(exonList) #length total
					infoTranscript[5]=len(exonList)-1#number junction CDS_CDS
					LenghtByTranscript[transcriptId+'-'+transcriptBiotype]=infoTranscript
			
	inputfile.close()
	return LenghtByTranscript
######################################################################################################################################################
def GetNumbersG4Transcript(filename, Coding):
	""" Create dictionary of number of G4 for each transcript depending the localisation
	if transcript is coding ,format =  [number G4 in section 5, cds, 3,exon, intron,junction 5, junction CDS, junction i, junction 3, total, junction cds_cds ]
	else format = [number G4 in exon, intron, junction CDS, junction i, total, junction cds_cds ]
	    Parameters
	    ----------
	     filename : string
		name of the file which contain information for each G4 in a transcript 
	    Coding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    Returns
	    -------
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	"""
	NumbersG4Transcript={}
	inputfile= open(filename,"r") # file opening for reading
	for line in inputfile: # for each transcript
		if ('InfoG4ByTranscript' not in line):
			words=line.split('\t')
			InfoG4ByTranscript=words[0].rstrip() 
			transcriptId=InfoG4ByTranscript.split('|')[0]
			chromosome=InfoG4ByTranscript.split('|')[1].split(':')[0]
			startG4=int(InfoG4ByTranscript.split('|')[1].split(':')[1].split('-')[0])
			endG4=int(InfoG4ByTranscript.split('|')[1].split(':')[1].split('-')[1])
			strand=InfoG4ByTranscript.split('|')[2]
			cGcC=words[1].rstrip() 
			G4Hunter=words[2].rstrip() 
			sequenceG4=words[3].rstrip() 
			G4NN=words[4].rstrip() 
			localisation=words[5].rstrip() 
			transcriptBiotype=words[6].rstrip() 
			key=transcriptId+'-'+transcriptBiotype
			if transcriptBiotype in Coding: # if the trancript is a protein coding 
				if (NumbersG4Transcript.has_key(key) == False): # if gene not contain G4 detected
						numberG4=[0,0,0,0,0,0,0,0,0,0,0] # number G4 in section 5, cds, 3,exon, intron,junction 5, junction CDS, junction i, junction 3, total, junction cds_cds
				else: # if gene contain G4 detected 
						numberG4=NumbersG4Transcript.get(key)
				
				position=GetPositionForLocalisation(localisation, 'Coding')
				numberG4[position]=numberG4[position]+1
				NumbersG4Transcript[key]=numberG4
				if (localisation != 'junction_CDS_CDS')	: # for total in transcript entire (all except junction_CDS_CDS) 
					numberG4[9]=numberG4[9]+1
					NumbersG4Transcript[key]=numberG4
				if (localisation == '5' or localisation == 'CDS' or localisation == '3'): #for total exon (5, CDS, 3)
					numberG4[3]=numberG4[3]+1
					NumbersG4Transcript[key]=numberG4
			else:
				if (NumbersG4Transcript.has_key(key) == False): # if gene not contain G4 detected
						numberG4=[0,0,0,0,0,0] # number G4 in exon, intron, junction CDS, junction i, total, junction cds_cds
				else: # if gene contain G4 detected 
						numberG4=NumbersG4Transcript.get(key)
				position=GetPositionForLocalisation(localisation,'LongNC' )
				numberG4[position]=numberG4[position]+1
				NumbersG4Transcript[key]=numberG4
				if (localisation != 'junction_ExonNC_ExonNC')	: # for total in transcript entire (all except junction_CDS_CDS) 
					numberG4[4]=numberG4[4]+1
					NumbersG4Transcript[key]=numberG4
	inputfile.close()
	return NumbersG4Transcript
######################################################################################################################################################
def GetFrequenceG4ByTranscript (NumbersG4Transcript,LenghtByTranscript, Coding):
	""" Create dictionary of frequence of G4 for each transcript by position
	if transcript is coding ,format frequence=  [frequence G4 in section 5, cds, 3,exon, intron,junction 5, junction CDS, junction i, junction 3, total, junction cds_cds]
	else format = [G4 in exon, intron, junction CDS, junction i, total, junction cds_cds]
	    Parameters
	    ----------
	    NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	    Coding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    Returns
	    ----------
	     FrequenceG4ByTranscript : dictionary
		 dictionary of frequence of G4 for each transcript by position
	"""
	FrequenceG4ByTranscript={}
	for key, numberG4 in NumbersG4Transcript.items(): # for each transcript
		
		lenghtG4=LenghtByTranscript.get(key) # ensembl of information by transcript , format [length 5,length cds,length 3,lenght exon,lenght intron,number junction 5_codant,number junction codant_Intron,number junction Intron_codant,number junction codant_3,length total, number junction CDS_CDS ]
		transcriptBiotype=key.split('-')[1]
		if (transcriptBiotype in Coding):
			frequencyG4=[0,0,0,0,0,0,0,0,0,0,0] # frequency by transcript
			positionConditionLenght=[0,1,2,3,4,9] #[length 5,length cds,length 3,lenght exon,lenght intron junction,length total ]
		else:
			frequencyG4=[0,0,0,0,0,0] # frequency by transcript , G4 in exon, intron, junction CDS, junction i, total, junction cds_cds
			positionConditionLenght=[0,1,4] #[lenght exon,lenght intron,length total ]
		for position in range(len(numberG4)):# for each length or portion (number junction)
			if (position in positionConditionLenght): # if lenght --> frequence per Kb
				multiplier=1000 # by kb
			else: # if number junction
				multiplier=1 #--> frequence par junction
			#~ print key, lenghtG4
			#~ print position, 'hello', lenghtG4[position]
			if (lenghtG4[position] != 0): # if the denominator is different to zero
				frequencyPosition=round((numberG4[position]/float(lenghtG4[position])*multiplier),6) # calcul the frequence
			else: # if the denominator is zero
				frequencyPosition=round(0,6) # frequence will be equal to 0
			frequencyG4[position]=frequencyPosition # place this frequence in list according to the position in the list
		FrequenceG4ByTranscript[key]=frequencyG4 # add the list of frequence in the dictionary (key = transcriptId+biotype of transcript)
	return FrequenceG4ByTranscript
######################################################################################################################################################
def GetTranscriptWithG4ByFamily(filename):
	""" Create dictionary , where for each family, is repertoried the transcript with contain a G4.
	    Parameters
	    ----------
	    filename : string
		name of the file which contain information for G4 in a transcript
	    Returns
	    ----------
	     TranscriptWithG4ByFamily : dictionary
		 dictionary of transcript which contain G4 for each family
	"""
	TranscriptWithG4ByFamily={}
	inputfile= open(filename,"r") # file opening for reading
	for line in inputfile: # for each transcript
		if ('InfoG4ByTranscript' not in line): # if line isn't the header
			words=line.split('\t')
			InfoG4ByTranscript=words[0].rstrip() 
			transcriptId=InfoG4ByTranscript.split('|')[0]
			transcriptBiotype=words[6].rstrip() 
			if (TranscriptWithG4ByFamily.has_key(transcriptBiotype) == False): 
				listeTranscript=[] 
			else: # if gene contain G4 detected 
				listeTranscript=TranscriptWithG4ByFamily.get(transcriptBiotype)
			if (transcriptId not in listeTranscript):
				listeTranscript.append(transcriptId)
				TranscriptWithG4ByFamily[transcriptBiotype]=listeTranscript	
	inputfile.close()
	return TranscriptWithG4ByFamily
######################################################################################################################################################
def GetTranscriptWithG4ByLocalisation(filename):
	""" Create dictionary , where for each family and localisation, is repertoried the transcript with contain a G4.
	    Parameters
	    ----------
	    filename : string
		name of the file which contain information for each G4 in a transcript
	    Returns
	    ----------
	     TranscriptWithG4ByLocalisation : dictionary
		 dictionary of transcript which contain G4 for each family and localisation
	"""
	TranscriptWithG4ByLocalisation={}
	inputfile= open(filename,"r") # file opening for reading
	for line in inputfile: # for each transcript
		if ('InfoG4ByTranscript' not in line): # if line isn't the header
			words=line.split('\t')
			InfoG4ByTranscript=words[0].rstrip() 
			transcriptId=InfoG4ByTranscript.split('|')[0]
			localisation=words[5].rstrip() 
			transcriptBiotype=words[6].rstrip() 
			if (TranscriptWithG4ByLocalisation.has_key(transcriptBiotype+'|'+localisation) == False):  # if first time we met the key
				listeTranscript=[] #list of transcript empty
			else: # if not first time 
				listeTranscript=TranscriptWithG4ByLocalisation.get(transcriptBiotype+'|'+localisation) # take the list already fill
			if (transcriptId not in listeTranscript): # if this transcript no contain in the list (no redondance needed)
				listeTranscript.append(transcriptId) # add the transcript to the list
				TranscriptWithG4ByLocalisation[transcriptBiotype+'|'+localisation]=listeTranscript # add the list in the dico	
			conditionExon=['5', 'CDS', '3','junction_5_CDS','junction_CDS_3'] # condition exonique
			if (localisation in conditionExon): # if the localisation is in the exon
				localisation='Exon' # change variable localisation
				if (TranscriptWithG4ByLocalisation.has_key(transcriptBiotype+'|'+localisation) == False): 
					listeTranscript=[] 
				else: 
					listeTranscript=TranscriptWithG4ByLocalisation.get(transcriptBiotype+'|'+localisation)
				if (transcriptId not in listeTranscript):
					listeTranscript.append(transcriptId)
					TranscriptWithG4ByLocalisation[transcriptBiotype+'|'+localisation]=listeTranscript # add info for exon
	inputfile.close()
	return TranscriptWithG4ByLocalisation
######################################################################################################################################################
def GetTranscriptByFamily(filename, BiotypeByTranscript, AnnotationTranscript):
	""" Create dictionary , where for each family, is repertoried the transcript (with or without G4)
	    Parameters
	    ----------
	    filename : string
		name of the file which contain information fr each transcript (except the transcript's biotype)
            BiotypeByTranscript : dico
		dictionary, where key =id of transcript and value = biotype of this transcript
	    AnnotationTranscript: dico
		dictionary with information of annotation for the transcripts By Ensembl (True if good anotation, else False)
	    Returns
	    ----------
	     TranscriptByFamily : dictionary
		 dictionary of transcript for  each family
	"""
	TranscriptByFamily={}
	inputfile= open(filename,"r") # file opening for reading
	for line in inputfile: # for each transcript
		words=line.split('|')
		transcriptId=words[0].rstrip() 
		transcriptBiotype=BiotypeByTranscript.get(transcriptId) # get biotype of this transcript
		annotationTranscript=AnnotationTranscript.get(transcriptId)
		if (annotationTranscript == True and transcriptBiotype is not None): # if annotation of this transcript is correct
			if (TranscriptByFamily.has_key(transcriptBiotype) == False): 
				listeTranscript=[] 
			else: 
				listeTranscript=TranscriptByFamily.get(transcriptBiotype)
			if (transcriptId not in listeTranscript):
				listeTranscript.append(transcriptId)
				TranscriptByFamily[transcriptBiotype]=listeTranscript	
	inputfile.close()
	return TranscriptByFamily
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


######################################################################################################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	parser.add_argument ('-p', '--path', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/mouse')
	parser.add_argument ('-CHR', '--CHROMOSOME', default = 'X')
	parser.add_argument ('-specie', '--specie', default = 'MM')
	parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-E', '--EXTENSION', default = 100)
	parser.add_argument ('-W', '--WINDOW', default = 60)
	parser.add_argument ('-S', '--STEP', default = 10)
	return parser
######################################################################################################################################################
def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path	# directory which contain all the directory CHROMOSOME
	CHROMOSOME=arg.CHROMOSOME	# CHROMOSOME to analyze
	specie=arg.specie	# specie to analyse
	THRESHOLD_G4H=arg.THRESHOLD_G4H	# threshold use to discriminate the score G4H (litterature = 0.9)
	THRESHOLD_CGCC=arg.THRESHOLD_CGCC	# threshold use to discriminate the score G4H (litterature = 4.5)
	THRESHOLD_G4NN=arg.THRESHOLD_G4NN	# threshold use to discriminate the score G4NN (litterature = 0.5)
	EXTENSION=arg.EXTENSION	# EXTENSION use for create the junction exon_exon previously
	WINDOW=arg.WINDOW	# size of window use in G4 screener
	STEP=arg.STEP	# size of step use in G4 screener


	#### variable for each family or ensembl
	Coding=['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene', 'IG_M_gene', 'IG_V_gene', 'IG_Z_gene', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay', 'protein_coding', 'TR_C_gene', 'TR_D_gene', 'TR_gene', 'TR_J_gene', 'TR_V_gene']
	Pseudogene=['transcribed_unitary_pseudogene','disrupted_domain', 'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene', 'IG_V_pseudogene', 'processed_pseudogene', 'pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene', 'translated_unprocessed_pseudogene', 'TR_J_pseudogene', 'TR_V_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene','polymorphic_pseudogene']
	LongNC=['macro_lncRNA','bidirectional_promoter_lncRNA','sense_intronic','3prime_overlapping_ncRNA','ambiguous_orf','antisense','lincRNA','ncrna_host','non_coding','processed_transcript','retained_intron','sense_overlapping']
	ShortNC=['scaRNA','miRNA', 'miRNA_pseudogene', 'misc_RNA', 'misc_RNA_pseudogene', 'Mt_rRNA' ,'Mt_tRNA', 'Mt_tRNA_pseudogene', 'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA', 'RNase_P_RNA', 'rRNA', 'rRNA_pseudogene', 'scRNA_pseudogene', 'snlRNA', 'snoRNA', 'snoRNA_pseudogene', 'snRNA', 'snRNA_pseudogene', 'SRP_RNA', 'tmRNA', 'tRNA', 'tRNA_pseudogene']
	Predictif=['TEC']
	NonCoding=['macro_lncRNA','bidirectional_promoter_lncRNA','sense_intronic','3prime_overlapping_ncRNA','ambiguous_orf','antisense','lincRNA','ncrna_host','non_coding','processed_transcript','retained_intron','sense_overlapping','3prime_overlapping_ncRNA','scaRNA','miRNA', 'miRNA_pseudogene', 'misc_RNA', 'misc_RNA_pseudogene', 'Mt_rRNA' ,'Mt_tRNA', 'Mt_tRNA_pseudogene', 'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA', 'RNase_P_RNA', 'rRNA', 'rRNA_pseudogene', 'scRNA_pseudogene', 'snlRNA', 'snoRNA', 'snoRNA_pseudogene', 'snRNA', 'snRNA_pseudogene', 'SRP_RNA', 'tmRNA', 'tRNA', 'tRNA_pseudogene']
	SupraFamily=[Coding, NonCoding, Pseudogene, Predictif]
	SupraFamilyName=['Coding', 'NonCoding', 'Pseudogene', 'Predictif']

	
	
	#### path for files depending of chromosome or for all
	if (CHROMOSOME == 'All'):	## if for all
		directory=path+'/All'
		index=directory+'/'+specie+'_transcript_unspliced_All_Index.txt'	
		indextranscriptBiotype=directory+'/'+specie+'_All_TranscriptType.txt'
		fileG4InTranscriptome=directory+'/'+specie+'_All_G4InTranscript.txt'
	else:	## if by CHROMOSOME
		directory=path+'/chr'+CHROMOSOME
		index=directory+'/'+specie+'_transcript_unspliced_chr'+CHROMOSOME+'_Index.txt'# file which contain info by transcript for this chromosome
		indextranscriptBiotype=path+'/transcriptType/transcriptType_chr'+CHROMOSOME+'.txt'	# file which contain biotype of transcript for this chromosome
		fileG4InTranscriptome=path+'/G4Screener/'+specie+'_chr'+CHROMOSOME+'_G4InTranscript.txt'
	

	#### Assignation of variable for the start
	BiotypeByTranscript={}	# dictionary of biotype for each transcript
	AnnotationTranscript={} # dictionary of annotation (True or False) for each transcript
	LenghtByTranscript={} # dictionary of lenght for each transcript
	NumbersG4Transcript={} 	 # dictionary of number of G4 for each transcript by position
	FrequenceG4ByTranscript={} # dictionary of frequence of G4 for each transcript by position
	TranscriptByFamily={} # dictionary of list of transcript for each biotype (Tr with or without G4)
	TranscriptWithG4ByFamily={} # dictionary of list of transcript for each biotype (only Tr with G4 )
	TranscriptWithG4ByLocalisation={} # dictionary of list of transcript which contain G4 for each biotype and localisation
	

	#### Filling dictionary
	BiotypeByTranscript=CreateDictionaryBiotypeByTranscript (indextranscriptBiotype) # Filling the dictionary of biotype for each transcript
	AnnotationTranscript=GetAnnotationTranscript(index,Coding,BiotypeByTranscript) # Filling the dictionary of annotation (True or False) for each transcript
	LenghtByTranscript=GetLenghtByTranscript(index, BiotypeByTranscript,AnnotationTranscript, Coding) # Filling the dictionary of lenght for each transcript
	NumbersG4Transcript=GetNumbersG4Transcript(fileG4InTranscriptome,Coding) # Filing the dictionary of number of G4 for each transcript by position
	FrequenceG4ByTranscript=GetFrequenceG4ByTranscript(NumbersG4Transcript,LenghtByTranscript, Coding)# Filing the dictionary of frequence of G4 for each transcript by position
	TranscriptByFamily=GetTranscriptByFamily(index, BiotypeByTranscript, AnnotationTranscript)# Filing the dictionary of list of transcript for each biotype (Tr with or without G4)
	TranscriptWithG4ByFamily=GetTranscriptWithG4ByFamily(fileG4InTranscriptome)# Filing the dictionary of list of transcript for each biotype (only Tr with G4 )
	TranscriptWithG4ByLocalisation=GetTranscriptWithG4ByLocalisation(fileG4InTranscriptome)# Filing the dictionary of list of transcript which contain G4 for each biotype and localisation
	




	path_out=directory+'/Subset' # directory pathway
	if os.path.isdir(path_out):## if directory output exist 	
		shutil.rmtree(path_out) # delate if --> refresh
	os.mkdir(path_out) # create directory output
	for family in SupraFamily:
		for subFamily in family:
			test=False
			if (TranscriptByFamily.has_key(subFamily) == True): 	# if this subset contain G4 detected
				if (subFamily in Coding):
					localisations=['5','3','CDS','Intron','junction_5_CDS','junction_CDS_3','junction_CDS_CDS','junction_CDS_Intron','junction_Intron_CDS'] # localisations possible for this subset
					typeFamily='Coding'
				else:
					localisations=['ExonNC','IntronNC','junction_ExonNC_ExonNC','junction_ExonNC_IntronNC','junction_IntronNC_ExonNC']
					typeFamily='Others'
				nbTRanscript=len(TranscriptByFamily[subFamily]) # nbr absoly of transcript for this subset
				for localisation in localisations: # for each localisation
					key=subFamily+'|'+localisation
					if (TranscriptWithG4ByLocalisation.has_key(key) == True and test==False):
						output= open(path_out+'/'+subFamily+'.csv',"w") # file opening for reading
						output.write('Localisation\tNbTRanscript\tNbTRanscriptWithG4\tPercent\tMean\tSD\tSE\tFrequency\n')
						test=True
				if (test==True):
					for localisation in localisations: # for each localisation
						key=subFamily+'|'+localisation # key for dictionary
						moy=var=sD=sE=0
						percent=0
						nbTRanscriptWithG4=0
						if (TranscriptWithG4ByLocalisation.has_key(key) == True): # if G4 detected in this family for this localisation
							nbTRanscriptWithG4=len(TranscriptWithG4ByLocalisation[key])# nbr absoly of transcript for this subset at this localisation
							percent=round(float(nbTRanscriptWithG4/float(nbTRanscript)*100),2) # percent of transcript with G4 at this localisation
							valuesWithG4=[] # list of value of frequence for this localisation and this subfamily, len(values) will be = to the nbTRanscriptWithG4
							position=GetPositionForLocalisation(localisation, typeFamily)# get position of this localisation
							for transcriptId in TranscriptWithG4ByLocalisation[key]: # for each transcript 
								valuesWithG4.append(FrequenceG4ByTranscript[transcriptId+'-'+subFamily][position]) # add the value for the frequence of this position
							valuesWithG4str = '\t'.join(str(e) for e in valuesWithG4)
							moy=round(mean(valuesWithG4),3) # mean of frequence
							var=CalculVariance(valuesWithG4, moy) # variance of mean
							sD=round(CalculSquareroot(var),3) # SD
							sE=round(CalculStandardError(sD, valuesWithG4),3) # SE
							output.write(localisation+'\t'+str(nbTRanscript)+'\t'+str(nbTRanscriptWithG4)+'\t'+str(percent)+'\t'+str(moy)+'\t'+str(sD)+'\t'+str(sE)+'\t'+str(valuesWithG4str)+'\n')
						else:
							output.write(localisation+'\t'+str(nbTRanscript)+'\t0\t0\t0\t0\t0\n')
				else:
					continue



	

	
	
	
	
	

main()
