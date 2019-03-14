#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
Copyright Université de Sherbrooke, Département de Biochimie et Département d'Informatique

contact: sarah.belhamiti@usherbrooke.ca

This software is a computer program whose calcul the statistic value of PG4r (density and percent)


---------------------------------------------------

``G4Calcul`` **module description**:

From a tabulate files which contain G4 regions detected in one chromosome (or all chromosomes), this module will study the the distribution of the regions in the transcriptome. 

.. moduleauthor:: Sarah.Belhamiti

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
		length=(abs(int(positionA)-int(positionB)))+1
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
			lengthCodant+=(GetLengthFraction(startExon,endExon))# add the lenght of this exon to the total (lengthCodant)
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
		elif(localisation == 'Exon'):
			position=3
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
		elif(localisation == 'Total'):
			position=9
		elif(localisation == 'junction_CDS_CDS'):
			position=10
		elif(localisation == 'Splicing'):
			position=11
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
		elif(localisation == 'Total'):
			position=4
		elif(localisation == 'junction_ExonNC_ExonNC'):
			position=5
		elif(localisation == 'Splicing'):
			position=6
		#### position 4 for total in transcript entire (all except junction_CDS_CDS))
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
		name of the file which contain informatiob fr each transcript (except the transcript's biotype)
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
					infoTranscript=[0,0,0,0,0,0,0,0,0,0,0,0] #format [length 5,length cds,length 3,lenght exon,lenght intron,number junction 5_codant,number junction codant_Intron,number junction Intron_codant,number junction codant_3,length total, number junction CDS_CDS ]
					infoTranscript[0]=GetLengthFraction(start5,end5) #length 5
					infoTranscript[1]=GetLengthCodant(exonList)-GetLengthFraction(start5,end5)-GetLengthFraction(start3,end3) #length cds
					infoTranscript[2]=GetLengthFraction(start3,end3) #length 3
					infoTranscript[3]=GetLengthCodant(exonList) #lenght exon
					if (intronList[0] != ''):
						infoTranscript[4]=GetLengthCodant(intronList) #lenght intron
					infoTranscript[5]=1 # number junction 5_codant
					infoTranscript[6]=len(exonList)-1 #  number junction codant_Intron
					infoTranscript[7]=len(exonList)-1 #number junction Intron_codant
					infoTranscript[8]=1 # number junction codant_3
					infoTranscript[9]=GetLengthTotal(exonList) #length total
					infoTranscript[10]=len(exonList)-1#number junction CDS_CDS
					infoTranscript[11]=infoTranscript[6]+infoTranscript[7]+infoTranscript[10] #(number junction on splicing events = CDS_Intron + INtronCDS+CDS_CDS)
					LenghtByTranscript[transcriptId+'-'+transcriptBiotype]=infoTranscript
					
				else:
					infoTranscript=[0,0,0,0,0,0,0] #format [lenght exon,lenght intron,number junction codant_Intron,number junction Intron_codant,length total, number junction CDS_CDS ]
					infoTranscript[0]=GetLengthCodant(exonList) #lenght exon
					if (intronList[0] != ''):
						infoTranscript[1]=GetLengthCodant(intronList) #lenght intron
					infoTranscript[2]=len(exonList)-1 #  number junction codant_Intron
					infoTranscript[3]=len(exonList)-1 #number junction Intron_codant
					infoTranscript[4]=GetLengthTotal(exonList) #length total
					infoTranscript[5]=len(exonList)-1#number junction CDS_CDS
					infoTranscript[6]=infoTranscript[2]+infoTranscript[3]+infoTranscript[5] #(number junction on splicing events = ExonNC_IntronNC + INtronNC_ExonNC +ExonNC_ExonNC)
					LenghtByTranscript[transcriptId+'-'+transcriptBiotype]=infoTranscript
			#if (len(exonList)==2):
				#print transcriptId,infoTranscript
			
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
						numberG4=[0,0,0,0,0,0,0,0,0,0,0,0] # number G4 in section 5, cds, 3,exon, intron,junction 5, junction CDS, junction i, junction 3, total, junction cds_cds, spicing
				else: # if gene contain G4 detected 
						numberG4=NumbersG4Transcript.get(key)
				
				position=GetPositionForLocalisation(localisation, 'Coding')
				numberG4[position]=numberG4[position]+1
				NumbersG4Transcript[key]=numberG4
				# for total in transcript entire (all) 
				numberG4[9]=numberG4[9]+1
				NumbersG4Transcript[key]=numberG4
				if (localisation == '5' or localisation == 'CDS' or localisation == '3' or localisation == 'junction_5_CDS' or localisation == 'junction_CDS_3'): #for total exon (5, CDS, 3, junction 5_CDS or  junction CDS_3)
					numberG4[3]=numberG4[3]+1
					NumbersG4Transcript[key]=numberG4
				if (localisation == 'junction_CDS_CDS' or localisation == 'junction_Intron_CDS' or localisation == 'junction_CDS_Intron'): #for total exon (5, CDS, 3, junction 5_CDS or  junction CDS_3)
					numberG4[11]=numberG4[11]+1
					NumbersG4Transcript[key]=numberG4
			else:
				if (NumbersG4Transcript.has_key(key) == False): # if gene not contain G4 detected
						numberG4=[0,0,0,0,0,0,0] # number G4 in exon, intron, junction CDS, junction i, total, junction cds_cds, splicing
				else: # if gene contain G4 detected 
						numberG4=NumbersG4Transcript.get(key)
				position=GetPositionForLocalisation(localisation,'LongNC' )
				numberG4[position]=numberG4[position]+1
				NumbersG4Transcript[key]=numberG4
				 # for total in transcript entire (all) 
				numberG4[4]=numberG4[4]+1
				NumbersG4Transcript[key]=numberG4
				if (localisation == 'junction_ExonNC_ExonNC' or localisation == 'junction_IntronNC_ExonNC' or localisation == 'junction_ExonNC_IntronNC'): #for total exon (5, CDS, 3, junction 5_CDS or  junction CDS_3)
					numberG4[6]=numberG4[6]+1
					NumbersG4Transcript[key]=numberG4
	inputfile.close()
	return NumbersG4Transcript
######################################################################################################################################################

def GetNumberTranscriptByCat (LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat,Coding):
	""" Create dictionary of number of transcript in each categorie 
	    Parameters
	    ----------
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 
	    Coding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    Returns
	    -------
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	"""
	for key, values in LenghtByTranscript.items():
		transcriptId=key.split('-')[0]
		transcriptBiotype=key.split('-')[1]
		if (NumberTranscriptByCat.has_key(transcriptBiotype) == False): # if biotype not contain in the dico
			if (transcriptBiotype in Coding):
				localisationName=['5','CDS','3','Exon','Intron','junction_5_CDS','junction_CDS_3','junction_CDS_Intron','junction_Intron_CDS','Total','junction_CDS_CDS','Splicing']
			else:
				localisationName=['ExonNC','IntronNC','junction_ExonNC_IntronNC','junction_IntronNC_ExonNC','Total','junction_ExonNC_ExonNC','Splicing']
			listeTranscript=[0]*(len(localisationName)+1)
		else: # if biotype already contained in the dico
			listeTranscript=NumberTranscriptByCat.get(transcriptBiotype)
		listeTranscript[0]=listeTranscript[0]+1 # add this transcipt in the number total of the biotype group
		if (NumbersG4Transcript.has_key(key) == True): # if this transcipt contain a G4 at an position
			numberG4=NumbersG4Transcript.get(key) # recup list of count G4 by position in the transcript
			for position in range(len(numberG4)): # for each position 
				if (numberG4[position] != 0): # if the transcript contain G4 at this precise position
					listeTranscript[position+1]= listeTranscript[position+1]+1 # add it at position +1
		NumberTranscriptByCat[transcriptBiotype]=listeTranscript
		
	return NumberTranscriptByCat
		
######################################################################################################################################################

def CreateRnaTypeDetected (LenghtByTranscript,rnaTypeDetected):
	""" Create list of RNA type where PG4r were detected 
	    Parameters
	    ----------
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     rnaTypeDetected : list
		list of RNA type where PG4r were detected
	    Returns
	    -------
	     rnaTypeDetected : list
		list of RNA type where PG4r were detected
	"""
	rnaTypeDetected=[]
	for key, values in LenghtByTranscript.items():
		biotypeTranscript=key.split('-')[1]
		if (biotypeTranscript not in rnaTypeDetected) :
			rnaTypeDetected.append(biotypeTranscript)
	return rnaTypeDetected

######################################################################################################################################################
def EnrichssmentByFamily (familyName, family,state,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat):
	""" Create list of calcul for a complete RNA family
	    Parameters
	    ----------
	    familyName : string
		Family where PG4r was detectd
	    family : list
		list of type contained in this family
	    state : string
		condirion where PG4r is detected (~ location)
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 

	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	if (familyName== 'Coding'):
		localisationName=['5','CDS','3','Exon','Intron','junction_5_CDS','junction_CDS_3','junction_CDS_Intron','junction_Intron_CDS','Total','junction_CDS_CDS']
	else:
		localisationName=['ExonNC','IntronNC','junction_ExonNC_IntronNC','junction_IntronNC_ExonNC','Total','junction_ExonNC_ExonNC']
	position=GetPositionForLocalisation(state, familyName)
	length=0
	num=0
	listeBiotype=[]
	for key, values in LenghtByTranscript.items(): # for each transcript in the transcriptom
		transcriptId=key.split('-')[0]
		transcriptBiotype=key.split('-')[1]
		if (transcriptBiotype in family): # if transcript in this family
			length=length+values[position]
			if (transcriptBiotype not in listeBiotype):
				listeBiotype.append(transcriptBiotype)
			if (NumbersG4Transcript.has_key(key)==True): # if transcript contain G4
				num=num+NumbersG4Transcript.get(key)[position]
	numTrAll=0
	numTrWithG4=0
	for biotype in listeBiotype:
		position=GetPositionForLocalisation(state, familyName)
		numTrAll=numTrAll+NumberTranscriptByCat.get(biotype)[0]
		numTrWithG4=numTrWithG4+NumberTranscriptByCat.get(biotype)[position+1]
	percent=round(numTrWithG4/float(numTrAll)*100,2)
	if ('junction' in state):
		return [familyName,length, num, round(num/float(length),4),numTrAll,numTrWithG4,percent]
	else:
		return [familyName,length, num, round(num/float(length)*1000,4),numTrAll,numTrWithG4,percent]
######################################################################################################################################################

def EnrichssmentByRnaType (familyName, family,state,LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat):
	""" Create list of calcul for a RNA type  
	    Parameters
	    ----------
	    familyName : string
		Family where PG4r was detectd
	    family : list
		list of type contained in this family
	    state : string
		condirion where PG4r is detected (~ location)
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	     rnaTypeDetected : list
		list of RNA type where PG4r were detected
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 

	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	if (familyName== 'Coding'):
		localisationName=['5','CDS','3','Exon','Intron','junction_5_CDS','junction_CDS_3','junction_CDS_Intron','junction_Intron_CDS','Total','junction_CDS_CDS']
	else:
		localisationName=['ExonNC','IntronNC','junction_ExonNC_IntronNC','junction_IntronNC_ExonNC','Total','junction_ExonNC_ExonNC']
	position=GetPositionForLocalisation(state, familyName)
	liste=[]
	for rnaType in family:
		if (rnaType in rnaTypeDetected):
			length=0
			num=0
			numTrAll=0
			numTrWithG4=0
			for key, values in LenghtByTranscript.items(): # for each transcript in the transcriptom
				transcriptId=key.split('-')[0]
				transcriptBiotype=key.split('-')[1]
				if (transcriptBiotype == rnaType): # if transcript in this family
					length=length+values[position]
					if (NumbersG4Transcript.has_key(key)==True): # if transcript contain G4
						num=num+NumbersG4Transcript.get(key)[position]
			position=GetPositionForLocalisation(state,familyName)
			numTrAll=NumberTranscriptByCat.get(rnaType)[0]
			numTrWithG4=NumberTranscriptByCat.get(rnaType)[position+1]
			if (length !=0):
				percent=round(numTrWithG4/float(numTrAll)*100,2)
				if ('junction' in state):
					liste.append([rnaType,length, num, round(num/float(length),4),numTrAll,numTrWithG4,percent])
				else:
					liste.append([rnaType,length, num, round(num/float(length)*1000,4),numTrAll,numTrWithG4,percent])
			else:
					liste.append([rnaType,length, num, round(num/float(length)*1000,4),numTrAll,numTrWithG4,0])
			
	return liste

######################################################################################################################################################


def EnrichssmentByRnaTypeLocalisation (familyName, family,state,LenghtByTranscript,NumbersG4Transcript, rnaType,NumberTranscriptByCat):	
	""" Create list of calcul for a RNA type  at a specific position 
	    Parameters
	    ----------
	    familyName : string
		Family where PG4r was detectd
	    family : list
		list of type contained in this family
	    state : string
		condirion where PG4r is detected (~ location)
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	    rnaTypeDetected : list
		list of RNA type where PG4r were detected
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 
	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	position=GetPositionForLocalisation(state, familyName)
	length=0
	num=0
	numTrAll=0
	numTrWithG4=0
	for key, values in LenghtByTranscript.items(): # for each transcript in the transcriptom
		transcriptId=key.split('-')[0]
		transcriptBiotype=key.split('-')[1]
		if (transcriptBiotype == rnaType): # if transcript in this family
			length=length+values[position]
			if (NumbersG4Transcript.has_key(key)==True): # if transcript contain G4
				num=num+NumbersG4Transcript.get(key)[position]

	position=GetPositionForLocalisation(state,familyName)
	numTrAll=NumberTranscriptByCat.get(rnaType)[0]
	numTrWithG4=NumberTranscriptByCat.get(rnaType)[position+1]

	percent=round(numTrWithG4/float(numTrAll)*100,2)

	if (length !=0):
		if ('junction' in state):
			return [rnaType,length, num, round(num/float(length),4),numTrAll,numTrWithG4,percent]
		else:
			return [rnaType,length, num, round(num/float(length)*1000,4),numTrAll,numTrWithG4,percent]
	else:
		return [rnaType,length, num, 0,numTrAll,numTrWithG4,percent]
######################################################################################################################################################

def EnrichssmentByFamilyLocalisation(familyName, family,state,LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat):
	""" Create ensembl of list of calcul for all RNA type  at a specific position for a category
	    Parameters
	    ----------
	    familyName : string
		Family where PG4r was detectd
	    family : list
		list of type contained in this family
	    state : string
		condirion where PG4r is detected (~ location)
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	    rnaTypeDetected : list
		list of RNA type where PG4r were detected
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 
	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	liste=[]
	position=GetPositionForLocalisation(state, familyName)
	for rnaType in family:
		if (rnaType in rnaTypeDetected):
			liste.append(EnrichssmentByRnaTypeLocalisation (familyName, family,state,LenghtByTranscript,NumbersG4Transcript, rnaType,NumberTranscriptByCat ))
	if ('junction' in state):
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste])),4)])
	else:
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste])),4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)]) # *100 because percent
	return liste
######################################################################################################################################################
def PrintInFile (outfilename,header,column,liste):
	""" print line for each element in a liste
	    Parameters
	    ----------
	    outfilename : string
		name file output
	    header : string
		header of the table
	    column : list
		list of column name
	    liste : list
		list of informations
	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	outfilename.write(header+'\n')
	outfilename.write(column+'\n')
	for element in liste:
		s=[str(item) for item in element]
		outfilename.write('\t'.join(s)+'\n')
	outfilename.write('\n')
				
######################################################################################################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	GITDIR=os.getcwd()+'/../'
	parser.add_argument ('-p', '--path', default = GITDIR+'results/')
	parser.add_argument ('-CHR', '--chromosome', default = 'all')
	parser.add_argument ('-specie', '--specie', default = 'HS')
	parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-E', '--EXTENSION', default = 100)
	parser.add_argument ('-W', '--WINDOW', default = 60)
	parser.add_argument ('-S', '--STEP', default = 10)
	parser.add_argument ('-c', '--choice', default = 'Coding')
	return parser
######################################################################################################################################################
def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path	# directory which contain all the directory chromosome
	chromosome=arg.chromosome	# chromosome to analyze
	specie=arg.specie	# specie to analyse
	THRESHOLD_G4H=arg.THRESHOLD_G4H	# threshold use to discriminate the score G4H (litterature = 0.9)
	THRESHOLD_CGCC=arg.THRESHOLD_CGCC	# threshold use to discriminate the score G4H (litterature = 4.5)
	THRESHOLD_G4NN=arg.THRESHOLD_G4NN	# threshold use to discriminate the score G4NN (litterature = 0.5)
	EXTENSION=arg.EXTENSION	# EXTENSION use for create the junction exon_exon previously
	WINDOW=arg.WINDOW	# size of window use in G4 screener
	STEP=arg.STEP	# size of step use in G4 screener
	choice=arg.choice


	GITDIR=os.getcwd()+'/../'

	#### variable for each family or ensembl
	Coding=['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene', 'IG_M_gene', 'IG_V_gene', 'IG_Z_gene', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay', 'protein_coding', 'TR_C_gene', 'TR_D_gene', 'TR_gene', 'TR_J_gene', 'TR_V_gene']
	Pseudogene=['transcribed_unitary_pseudogene','disrupted_domain', 'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene', 'IG_V_pseudogene', 'processed_pseudogene', 'pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene', 'translated_unprocessed_pseudogene', 'TR_J_pseudogene', 'TR_V_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene','polymorphic_pseudogene']
	LongNC=['macro_lncRNA','bidirectional_promoter_lncRNA','sense_intronic','3prime_overlapping_ncRNA','ambiguous_orf','antisense','lincRNA','ncrna_host','non_coding','processed_transcript','retained_intron','sense_overlapping']
	ShortNC=['vaultRNA','scaRNA','miRNA', 'miRNA_pseudogene', 'misc_RNA', 'misc_RNA_pseudogene', 'Mt_rRNA' ,'Mt_tRNA', 'Mt_tRNA_pseudogene', 'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA', 'RNase_P_RNA', 'rRNA', 'rRNA_pseudogene','scRNA', 'scRNA_pseudogene', 'snlRNA', 'snoRNA', 'snoRNA_pseudogene', 'snRNA', 'snRNA_pseudogene', 'SRP_RNA', 'tmRNA', 'tRNA', 'tRNA_pseudogene','ribozyme']
	Predictif=['TEC']
	NonCoding=['vaultRNA','macro_lncRNA','bidirectional_promoter_lncRNA','sense_intronic','3prime_overlapping_ncRNA','ambiguous_orf','antisense','lincRNA','ncrna_host','non_coding','processed_transcript','retained_intron','sense_overlapping','3prime_overlapping_ncRNA','scaRNA','miRNA', 'miRNA_pseudogene', 'misc_RNA', 'misc_RNA_pseudogene', 'Mt_rRNA' ,'Mt_tRNA', 'Mt_tRNA_pseudogene', 'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA', 'RNase_P_RNA', 'rRNA', 'rRNA_pseudogene','scRNA','scRNA_pseudogene', 'snlRNA', 'snoRNA', 'snoRNA_pseudogene', 'snRNA', 'snRNA_pseudogene', 'SRP_RNA', 'tmRNA', 'tRNA', 'tRNA_pseudogene','ribozyme']
	SupraFamily=[Coding, NonCoding, Pseudogene, Predictif]
	SupraFamilyName=['Coding', 'NonCoding', 'Pseudogene', 'Predictif']

	
	
	#### path for files depending of chromosome or for all
	if (chromosome == 'all'):	## if for all
		directory=path+'/all'
		index=directory+'/'+specie+'_transcript_unspliced_All_Index.txt'	
		indextranscriptBiotype=directory+'/transcriptType_All'
		fileG4InTranscriptome=directory+'/HS_All_G4InTranscript.txt'
	else:	## if by chromosome
		directory=path+'/chr'+chromosome
		index=directory+'/'+specie+'_transcript_unspliced_chr'+chromosome+'_Index.txt'# file which contain info by transcript for this chromosome
		indextranscriptBiotype=path+'/transcriptType/transcriptType_chr'+chromosome	# file which contain biotype of transcript for this chromosome
		fileG4InTranscriptome='/home/local/USHERBROOKE/bels2814/Documents/These/HUMAN/HS_chr'+chromosome+'_G4InTranscript.txt'



	#### Assignation of variable for the start
	BiotypeByTranscript={}	# dictionary of biotype for each transcript
	AnnotationTranscript={} # dictionary of annotation (True or False) for each transcript
	LenghtByTranscript={} # dictionary of lenght for each transcript
	NumbersG4Transcript={} 	 # dictionary of number of G4 for each transcript by position
	NumberTranscriptByCat={}# dictionary of list of number of each for each rna type (position 0 = all gene, position 1 = gene with G4)
	rnaTypeDetected=[] # list of rnaType detected in this species
	

	#### Filling dictionary
	BiotypeByTranscript=CreateDictionaryBiotypeByTranscript (indextranscriptBiotype) # Filling the dictionary of biotype for each transcript
	AnnotationTranscript=GetAnnotationTranscript(index,Coding,BiotypeByTranscript) # Filling the dictionary of annotation (True or False) for each transcript
	LenghtByTranscript=GetLenghtByTranscript(index, BiotypeByTranscript,AnnotationTranscript, Coding) # Filling the dictionary of lenght for each transcript
	NumbersG4Transcript=GetNumbersG4Transcript(fileG4InTranscriptome,Coding) # Filing the dictionary of number of G4 for each transcript by position
	NumberTranscriptByCat=GetNumberTranscriptByCat(LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat,Coding)
	



	
	rnaTypeDetected=CreateRnaTypeDetected (LenghtByTranscript,rnaTypeDetected)
	
	column='RNA Type\tTotal Length\t#PG4r\tDensity\t#Tr\t#TrWithPG4r\t%Tr'
	outfilename= open(GITDIR+'resultsTable/'+choice+'.csv',"w")
	
	if (choice == 'Coding'):
		liste=[]
		condition='Total'
		liste.append(EnrichssmentByFamily ('Coding', Coding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('NonCoding', NonCoding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Pseudogene', Pseudogene,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Predictif', Predictif,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header='Whole transcriptome'
		PrintInFile (outfilename,header,column,liste)
		###
		liste=EnrichssmentByRnaType ('Coding', Coding,'Total',LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header='By mRNA Subclass '
		PrintInFile (outfilename,header,column,liste)
		###
		conditions=[['Exon','Intron','Splicing'],['5','CDS','3', 'junction_5_CDS', 'junction_CDS_3'],['junction_CDS_Intron','junction_Intron_CDS','junction_CDS_CDS']]
		headers=[['Exon','Intron','Splicing Site'],["5'UTR",'CDS',"3'UTR", 'Codon start', 'Codon stop'],['Donor splicing site','Acceptor splicing site','Splicing junction']]
		for condition in conditions:
			for state in condition:
				#liste=EnrichssmentByFamily ('Coding', Coding,state,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat)
				liste=EnrichssmentByFamilyLocalisation('Coding', Coding,state,LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
				header=headers[conditions.index(condition)][condition.index(state)]
				PrintInFile (outfilename,header,column,liste)

	elif (choice == 'NonCoding'):
		liste=[]
		condition='Total'
		liste.append(EnrichssmentByFamily ('Coding', Coding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('NonCoding', NonCoding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Pseudogene', Pseudogene,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Predictif', Predictif,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header='Whole transcriptome'
		PrintInFile (outfilename,header,column,liste)
		###
		liste.append(EnrichssmentByFamily ('LongNC', LongNC,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('ShortNC', ShortNC,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header='By ncRNA groupe'
		PrintInFile (outfilename,header,column,liste)
		###
		liste=EnrichssmentByRnaType ('LongNC', LongNC,'Total',LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header=header='By long ncRNA Subclass '
		PrintInFile (outfilename,header,column,liste)
		###
		conditions=[['ExonNC','IntronNC','Splicing'],['junction_ExonNC_IntronNC','junction_IntronNC_ExonNC','junction_ExonNC_ExonNC']]
		headers=[['Exon','Intron','Splicing Site'],['Donor splice site','Acceptor splice site','Splicing Junction ']]
		for condition in conditions:
			for state in condition:
				#liste=EnrichssmentByFamily ('LongNC', LongNC,state,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat)
				liste=EnrichssmentByFamilyLocalisation('LongNC', LongNC,state,LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
				header=headers[conditions.index(condition)][condition.index(state)]
				PrintInFile (outfilename,header,column,liste)
		###
		liste=EnrichssmentByRnaType ('ShortNC', ShortNC,'Total',LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header=header='By short ncRNA Subclass '
		PrintInFile (outfilename,header,column,liste)
	
	
	else:
		

		liste=[]
		condition='Total'
		liste.append(EnrichssmentByFamily ('Coding', Coding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('NonCoding', NonCoding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Pseudogene', Pseudogene,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Predictif', Predictif,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])

		header='Whole transcriptome'
		PrintInFile (outfilename,header,column,liste)
		
			
		#print '-----------------------------------------------------'


		liste=EnrichssmentByRnaType ('Pseudogene', Pseudogene,'Total',LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])



		header='By pseudogene Subclass '
		PrintInFile (outfilename,header,column,liste)
		####
		conditions=[['ExonNC','IntronNC','Splicing'],['junction_ExonNC_IntronNC','junction_IntronNC_ExonNC','junction_ExonNC_ExonNC']]
		headers=[['Exon','Intron','Splicing Site'],['Donor splice site','Acceptor splice site','Splicing Junction ']]
		for condition in conditions:
			
			for state in condition:
				#liste=EnrichssmentByFamily ('Pseudogene', Pseudogene,state,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat)
				
				liste=EnrichssmentByFamilyLocalisation('Pseudogene', Pseudogene,state,LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
				header=headers[conditions.index(condition)][condition.index(state)]
				PrintInFile (outfilename,header,column,liste)
				
		
	
	


main()
