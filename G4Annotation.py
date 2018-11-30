#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
Copyright Université de Sherbrooke, Département de Biochimie et Département d'Informatique

contact: sarah.belhamiti@usherbrooke.ca

This software is a computer program whose annote the G4 region of one chromosome of one specie.


---------------------------------------------------

``G4Annotation`` **module description**:

From a csv files which contain windows of G4 screener of a chromosome, this module annote
the G4 region for each transcript. The program create differents files output:

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

######################################################################################################################################################
def mean(liste): # give mean of a list
	""" Give mean of a list
	    Parameters
	    ----------
	    liste : list
		list of values
	    Returns
	    -------
	    sum(liste)/len(liste)
		float
		mean of a list
	"""
   	return sum(liste)/len(liste)

######################################################################################################################################################
def positionChromosomiqueGenePositif(position, EXTENSION, startIntron, endIntron): # convert a position of an Intron in chromosomique position (gene +)
	""" Give real position (chromosomique) of start or end of G4 around a junction from coordination of an intron FOR A GENE POSITIF
	    Parameters
	    ----------
	    position : integer
		position of startG4 or end G4 give by G4 screener, not about the chromosome
            EXTENSION : integer constante
		extension used for create sequence of junctions from coordination of each intron
            startIntron : integer
		start of intron, used for this junction
	    endIntron : integer
		end of intron, used for this junction
		
	    Returns
	    -------
	    position
		integer
		position chromosomique of the start or end G4
	"""
	if (int(position) <= EXTENSION):	## because start G4 classifier from 0
		position=int(startIntron)-1-EXTENSION+int(position)
	else:	# if from sequence aval
		position=int(endIntron)+1-EXTENSION+int(position)-1
	return position
######################################################################################################################################################
def positionChromosomiqueGeneNegatif(position, EXTENSION, startIntron, endIntron):	# convert a position of an Intron in chromosomique position (gene -)
	""" Give real position (chromosomique) of start or end of G4 around a junction from coordination of an intron FOR A GENE NEGATIF
	    Parameters
	    ----------
	    position : integer
		position of startG4 or end G4 give by G4 screener, not about the chromosome
            EXTENSION : integer constante
		extension used for create sequence of junctions from coordination of each intron
            startIntron : integer
		start of intron, used for this junction
	    endIntron : integer
		end of intron, used for this junction
		
	    Returns
	    -------
	    position
		integer
		position chromosomique of the start or end G4
	"""
	if (int(position) <= EXTENSION):	## because strart G4 classifier from 0
		position=int(startIntron)+1+EXTENSION-int(position)
	else:	# if from sequence aval
		position=int(endIntron)+EXTENSION-int(position)
	return position
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
		transcript=words[1].rstrip()	#ID of transcript
		biotypeTranscript=words[3].rstrip() # biotype of the transcript
		if (dico.has_key(transcript) == False):	# if transcript not contain his biotype in the dico
			dico[transcript]=biotypeTranscript # create an entry (biotype) for this transcript
	return dico	# return dico with all transcript of this chromosome
######################################################################################################################################################
def CreateDictionaryStrandByGene (filename): 
	""" Create dictionary with as value the strand for all genes of the chromosome

	    Parameters
	    ----------
	    filename : string
		name of file which contain the liste of informations for each transcript of this chromosome

	    Returns
	    -------
	    dico
		dictionary
		where key = ID of gene, value = strand of this gene
	"""
	dico={}	#Creation dictionary empty, where will be : key = gene ID, value = strand of this gene
	inputfile= open(filename,"r")	# file opening for reading
	for line in inputfile:	# for each line in the file (for each transcript)
		words=line.split('|')	# parsing by the separator, here '|'
		gene=words[1].rstrip()	# gene ID
		strand=words[3].rstrip() # strand of this gene
		if (dico.has_key(gene) == False):	# if gene not contain his strand in the dico
			dico[gene]=strand # create an entry (strand) for this gene
	return dico	# return dico with all gene of this chromosome

######################################################################################################################################################
def ReturnG4InGene(G4DetectedInGene, inputfile, parametersTool, StrandByGene): ## for G4 in gene
	""" Add informations in dictionary with informations of region G4 detected (score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		regions of G4 discovered in the genes

	    Parameters
	    ----------
	    G4DetectedInGene : dictionary 
		dictionary will contain the information of G4 (score cgCc, scote G4Hunter, sequence and score G4NN) for all the region G4 detected
	    inputfile : string
		name of file which contain the liste of informations for each transcript of this chromosome
	    parametersTool : list of value use by G4 screener
		parametersTool[0]=float constante,threshold use to discriminate the score cGcC
		parametersTool[1]= float constante,threshold use to discriminate the score G4H
		parametersTool[2]= float constante,threshold use to discriminate the score G4NN
		parametersTool[3]= integer constante,size of window use in G4 screener	 
		parametersTool[4]=integer constante,size of step use in G4 screener
	    StrandByGene : dictionary
		dictionary with as value the strand for all genes of the chromosome

	    Returns
	    -------
	    G4DetectedInGene
		dictionary with informations of region G4 detected (score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		regions of G4 discovered in the genes
	"""
	THRESHOLD_CGCC=parametersTool[0]
	THRESHOLD_G4H=parametersTool[1]
	THRESHOLD_G4NN=parametersTool[2]
	WINDOW=parametersTool[3]
	STEP=parametersTool[4]
	oldPassed=False
	passed= False
	inputfile= open(inputfile,"r")	# file opening for reading
	for line in inputfile:	# for each file csv create by G4 screener
		if (re.search('^[0-9]', line)): # if the line is not the header of the file 
			words=line.split('\t') # parsing by the separator, here '\t' (file csv)
			numRow=words[0].rstrip() #numero of row for this for this window
			description=words[1].rstrip() # description of the sequence ( geneId|startBorder|endBorder )
			gene=words[1].rstrip().split("|")[0] # geneID
			strand=StrandByGene.get(gene) # get strand for this gene from the dictionary StrandByGene
			startBorder=int(words[1].rstrip().split("|")[1]) # start of the gene sequence
			endBorder=int(words[1].rstrip().split("|")[2])  # end of the gene sequence
			cGcC=float(words[2].rstrip()) #  score cGcC for this window
			g4H=float(words[3].rstrip()) #  score G4Hunter for this window 
			sequence=words[4].rstrip()	# sequence for this window
			startWindow=int(words[5].rstrip())	# start window from G4screener
			endWindow=int(words[6].rstrip()) ## end window from G4screener
			g4NN=float(words[7].rstrip()) # score G4NN for this window
			if (cGcC>=THRESHOLD_CGCC and g4H>=THRESHOLD_G4H and g4NN>=THRESHOLD_G4NN): # if window contain value > several threshold
				passed=True # passage over the differents threshold for this window
				if (oldPassed != passed): # if it's the first windows (beginning if the passage)
					listeCGcC=[] # list which will contain the values of score cgCG for windows over the THRESHOLD_CGCC
					listeG4H=[]  # list which will contain the values of score G4H for windows over the THRESHOLD_G4H
					listeG4NN=[]  # list which will contain the values of score G4NN for windows over the THRESHOLD_G4NN
					descriptionOverThreshold=description # assignation of description for this variable
					gene_startG4=gene
					sequenceG4=sequence
					oldPassed=passed # assignation, we pass over the thresholds
					if (strand == str(1)): # if gene positif
						startG4=startWindow # startG4 same as startWindow of G4 screener
						endG4=endWindow # endG4 same as endWindow of G4 screener
					else:	# if gene negatif
						startG4=int(endBorder)-(int(startWindow)-int(startBorder)) # calcul of startG4
						endG4=int(startG4)-(len(sequence))+1 # calcul of endG4
					listeCGcC=[cGcC] # add score cGCc  in the list 
					listeG4Hunter=[g4H] # add score G4H  in the list 
					listeG4NN=[g4NN] # add score G4NN  in the list 	
				else: # if it's not the first windows above the thresholds
					if (len(sequence) < WINDOW): 
					# if the sequence is shorter than the sequence of the windows (generraly last windows above the thresolds)
						sequenceG4=sequenceG4+sequence[-(len(sequence)-(WINDOW-STEP)):] # if last window for gene
					else: #
						sequenceG4=sequenceG4+sequence[-STEP:] # take the stepsize added

					if (strand == str(1)): # if gene positif
						endG4=endWindow # endG4 moves from the end of the first window,  same as endWindow of G4 screener
					else: # if gene negatif
						endG4=int(startG4)-(len(sequenceG4))+1 # endG4 moves from the end of the first window, # calcul of endG4
					lastRow=numRow # numRow of the last windows study over the threshols
					listeCGcC.append(cGcC) # add score cGcC for the windows over the threshold
					listeG4Hunter.append(g4H) # add score G4H for the windows over the threshold
					listeG4NN.append(g4NN) # add score G4NN for the windows over the threshold
			if (cGcC<THRESHOLD_CGCC or g4H<THRESHOLD_G4H or g4NN<THRESHOLD_G4NN or descriptionOverThreshold != description):
			# if one of the score is under his threshold or if this windows contain info of an other gene
				passed=False # passed passed became false (pass under the thresolds)
				if (oldPassed != passed ): # last windows before under the thresolds
					meanCGcC=mean(listeCGcC) # mean of the score cGcC will be the score for this region containing G4
					meanG4Hunter=mean(listeG4Hunter) # mean of the score G4Hunter will be the score for this region containing G4
					meanG4NN=mean(listeG4NN) # mean of the score G4NN will be the score for this region containing G4
					oldPassed=passed # assignation, we pass under the thresholds
					if (strand == str(1)):
						headerG4=gene+"|"+str(startG4)+"|"+str(endG4)+"|"+strand
					elif (strand == str(-1)):
						headerG4=gene+"|"+str(endG4)+"|"+str(startG4)+"|"+strand
					else:
						continue #strand == None, some gene doesnt have annotation
					if (G4DetectedInGene.has_key(headerG4) == False and strand != None): 
						G4DetectedInGene[headerG4]=str(meanCGcC), str(meanG4Hunter),sequenceG4 , str(meanG4NN)	
	return G4DetectedInGene
######################################################################################################################################################
def G4IsOnJunction(startG4, endG4, startBorder, endBorder):
	onJunction=False
	if ((startG4 <startBorder and endG4 >startBorder) or (startG4 >startBorder and endG4<startBorder)):
		onJunction=True
		
	return onJunction
######################################################################################################################################################
def ReturnG4InJunction(G4DetectedInJunction, inputfile, parametersTool,EXTENSION, StrandByGene ): ## for G4 in junction 
	""" Add informations in dictionary with informations of region G4 detected (score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		junctions of G4 discovered in the genes

	    Parameters
	    ----------
	    G4DetectedInJunction : dictionary 
		dictionary will contain the information of G4 (score cgCc, scote G4Hunter, sequence and score G4NN) for all the region G4 detected
	    inputfile : string
		name of file which contain the liste of informations for each transcript of this chromosome
	    parametersTool : list of value use by G4 screener
		parametersTool[0]=float constante,threshold use to discriminate the score cGcC
		parametersTool[1]= float constante,threshold use to discriminate the score G4H
		parametersTool[2]= float constante,threshold use to discriminate the score G4NN
		parametersTool[3]= integer constante,size of window use in G4 screener	 
		parametersTool[4]=integer constante,size of step use in G4 screener
	    EXTENSION : integer constante
		extension use for create the junction exon_exon or csd_cds previously
	    StrandByGene : dictionary
		dictionary with as value the strand for all genes of the chromosome

	    Returns
	    -------
	    G4DetectedInJunction
		dictionary with informations of region G4 detected (score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		junctions of G4 discovered in the genes
	"""
	THRESHOLD_CGCC=parametersTool[0] #
	THRESHOLD_G4H=parametersTool[1]
	THRESHOLD_G4NN=parametersTool[2]
	WINDOW=parametersTool[3]
	STEP=parametersTool[4]
	oldPassed=False
	passed= False
	descriptionOverThreshold=''
	inputfile= open(inputfile,"r")	# file opening for reading
	for line in inputfile:	# for each file csv create by G4 screener
		if (re.search('^[0-9]', line)): # if the line is not the header of the file 
			words=line.split('\t') # parsing by the separator, here '\t' (file csv)
			numRow=words[0].rstrip() #numero of row for this for this window
			description=words[1].rstrip() # description of the sequence ( geneId|startBorder|endBorder )
			gene=words[1].rstrip().split("|")[0] # geneID
			strand=StrandByGene.get(gene) # get strand for this gene from the dictionary StrandByGene
			startBorder=int(words[1].rstrip().split("|")[1].split("-")[0]) # start of the gene sequence
			endBorder=int(words[1].rstrip().split("|")[1].split("-")[1])  # end of the gene sequence
			cGcC=float(words[2].rstrip()) #  score cGcC for this window
			g4H=float(words[3].rstrip()) #  score G4Hunter for this window 
			sequence=words[4].rstrip()	# sequence for this window
			startWindow=int(words[5].rstrip())	# start window from G4screener
			endWindow=int(words[6].rstrip()) ## end window from G4screener
			g4NN=float(words[7].rstrip()) # score G4NN for this window
			if (cGcC>=THRESHOLD_CGCC and g4H>=THRESHOLD_G4H and g4NN>=THRESHOLD_G4NN and startWindow): # if window contain value > several threshold
				onJunction=False
				passed=True # passage over the differents threshold for this window
				if (oldPassed != passed): # if it's the first windows (beginning if the passage)
					listeCGcC=[] # list which will contain the values of score cgCG for windows over the THRESHOLD_CGCC
					listeG4H=[]  # list which will contain the values of score G4H for windows over the THRESHOLD_G4H
					listeG4NN=[]  # list which will contain the values of score G4NN for windows over the THRESHOLD_G4NN
					descriptionOverThreshold=description # assignation of description for this variable
					gene_startG4=gene
					sequenceG4=sequence
					startFirstWindow=startBorder
					endFirstWindow=endBorder
					oldPassed=passed # assignation, we pass over the thresholds
					#if (strand == str(1)): # if gene positif
					startG4=int(startWindow) # startG4 same as startWindow of G4 screener
					endG4=int(endWindow) # endG4 same as endWindow of G4 screener
					#else:	# if gene negatif
						#startG4=int(endBorder)-(int(startWindow)-int(startBorder)) # calcul of startG4
						#endG4=int(startG4)-(len(sequence))+1 # calcul of endG4
					listeCGcC=[cGcC] # add score cGCc  in the list 
					listeG4Hunter=[g4H] # add score G4H  in the list 
					listeG4NN=[g4NN] # add score G4NN  in the list 	
				else: # if it's not the first windows above the thresholds
					if (len(sequence) < WINDOW): 
					# if the sequence is shorter than the sequence of the windows (generraly last windows above the thresolds)
						sequenceG4=sequenceG4+sequence[-(len(sequence)-(WINDOW-STEP)):] # if last window for gene
					else: #
						sequenceG4=sequenceG4+sequence[-STEP:] # take the stepsize added
					#if (strand == str(1)): # if gene positif
					endG4=endWindow # endG4 moves from the end of the first window,  same as endWindow of G4 screener
					#else: # if gene negatif
						#endG4=int(startG4)-(len(sequenceG4))+1 # endG4 moves from the end of the first window, # calcul of endG4
					lastRow=numRow # numRow of the last windows study over the threshols
					listeCGcC.append(cGcC) # add score cGcC for the windows over the threshold
					listeG4Hunter.append(g4H) # add score G4H for the windows over the threshold
					listeG4NN.append(g4NN) # add score G4NN for the windows over the threshold
				
			if (cGcC<THRESHOLD_CGCC or g4H<THRESHOLD_G4H or g4NN<THRESHOLD_G4NN or descriptionOverThreshold != description):
			# if one of the score is under his threshold or if this windows contain info of an other gene
				passed=False # passed passed became false (pass under the thresolds)
				if (oldPassed != passed ): # last windows before under the thresolds
					meanCGcC=mean(listeCGcC) # mean of the score cGcC will be the score for this region containing G4
					meanG4Hunter=mean(listeG4Hunter) # mean of the score G4Hunter will be the score for this region containing G4
					meanG4NN=mean(listeG4NN) # mean of the score G4NN will be the score for this region containing G4
					oldPassed=passed # assignation, we pass under the thresholds
					if (strand == str(1)): # if gene positif
						startG4=positionChromosomiqueGenePositif(startG4, EXTENSION, startFirstWindow, endFirstWindow)	
						endG4=positionChromosomiqueGenePositif(endG4, EXTENSION, startFirstWindow, endFirstWindow)
					else: # if gene negatif
						startG4=positionChromosomiqueGeneNegatif(startG4, EXTENSION, startFirstWindow, endFirstWindow)	
						endG4=positionChromosomiqueGeneNegatif(endG4, EXTENSION, startFirstWindow, endFirstWindow)
					if (strand == str(1)):
						headerG4=gene+"|"+str(startG4)+"|"+str(endG4)+"|"+strand
					elif (strand == str(-1)):
						headerG4=gene+"|"+str(endG4)+"|"+str(startG4)+"|"+strand
					else:
						continue #strand == None, some gene doesnt have annotation
					onJonction=G4IsOnJunction(startG4, endG4, startBorder, endBorder)	
					if (G4DetectedInJunction.has_key(headerG4) == False and onJonction==True):
						G4DetectedInJunction[headerG4]=str(meanCGcC), str(meanG4Hunter),sequenceG4 , str(meanG4NN)
					
	return G4DetectedInJunction	

######################################################################################################################################################
def BorderOfTranscript(start5 , end5 , start3 , end3, exonList, intronList, strand):
	""" Give coordinate of borde for a transcript
	    Parameters
	    ----------
	    start5 : string 
		coordinate of beginning for the region 5', empty is transcript doesn't have region 5'
	    end5 : string
		coordinate of ending for the region 5', empty is transcript doesn't have region 5'
	    start3 : string 
		coordinate of beginning for the region 3', empty is transcript doesn't have region 3'
	    end3 : string
		coordinate of ending for the region 3', empty is transcript doesn't have region 3'
	    exonList : list 
 		list of exons contained in this transcript
	    intronList : list
		list of introns contained in this transcript
	    strand : string
		representation of the gene sens (1 for gene positif, -1 for gene negatif)
	    Returns
	    -------
	    [borderInferior,borderSuperior]
		liste of couple coordinate of border for a transcript
		borderInferior== border near to 5'
		borderSuperior== border near to 3'
	"""
	listBorder=[]
	if (start5 == ""):	# if transcript doesn't have 5'
		if (end3 == ""): # if transcript doesn't have 3'
			# border will be given by exons and introns
			for couple_exon in exonList:	# for each couple of exon of this transcript
				if not (couple_exon==""):
					startExon=couple_exon.split("-")[0] # start of this exon
					endExon=couple_exon.split("-")[1] # start of this exon
					listBorder.append(int(startExon))	# add information in list
					listBorder.append(int(endExon))	# add information in list
			for couple_Intron in intronList: # for each couple of intron of this transcript
				if not (couple_Intron==""):
					startIntron=couple_Intron.split("-")[0]	# start of this intron
					endIntron=couple_Intron.split("-")[1]	# end of this intron
					listBorder.append(int(startIntron))	# add information in list
					listBorder.append(int(endIntron))	# add information in list
			if (strand == str(1)):	#if gene positif, 5'< 3'
				borderInferior=min(listBorder)
				borderSuperior=max(listBorder)
			else: 	#if gene negatif, 5'> 3'
				borderInferior=max(listBorder)
				borderSuperior=min(listBorder)
		else: # if transcript doesn't have 5'
			# border will be given by exons and introns, and 3'
			for couple_exon in exonList:	# for each couple of exon of this transcript
				if not (couple_exon==""):
					startExon=couple_exon.split("-")[0] # start of this exon
					endExon=couple_exon.split("-")[1] # start of this exon
					listBorder.append(int(startExon))	# add information in list
					listBorder.append(int(endExon))	# add information in list
			for couple_Intron in intronList: # for each couple of intron of this transcript
				if not (couple_Intron==""):
					startIntron=couple_Intron.split("-")[0]	# start of this intron
					endIntron=couple_Intron.split("-")[1]	# end of this intron
					listBorder.append(int(startIntron))	# add information in list
					listBorder.append(int(endIntron))	# add information in list
			if (strand == str(1)): # if gene positif
				borderInferior=min(listBorder)
			else:	# if gene negatif
				borderInferior=max(listBorder)
			borderSuperior=int(end3)
	elif (start3 == ""): ##  if transcript doesn't have 3'
		# border will be given by exons and introns and 5'
		listBorder.append(start5)
		listBorder.append(end5)
		for couple_exon in exonList:
			if not (couple_exon==""):
				startExon=couple_exon.split("-")[0]
				endExon=couple_exon.split("-")[1]
				listBorder.append(int(startExon))
				listBorder.append(int(endExon))
		for couple_Intron in intronList:
			if not (couple_Intron==""):
				startIntron=couple_Intron.split("-")[0]
				endIntron=couple_Intron.split("-")[1]
				listBorder.append(int(startIntron))
				listBorder.append(int(endIntron))
		if (strand == str(1)):
			borderSuperior=max(listBorder)
		else:
			borderSuperior=min(listBorder)
		borderInferior=int(start5)
	else: ##  if transcript has 5', 3', exons and introns
		borderInferior=int(start5)
		borderSuperior=int(end3)
	return [borderInferior,borderSuperior]
######################################################################################################################################################
def G4IsInTranscript(strand, coordG4, borderTranscript):	
	""" Give boolean, if a G4 is located in a transcript
	    Parameters
	    ----------
	    strand : string
		representation of the gene sens (1 for gene positif, -1 for gene negatif)
	    coordG4 : couple of integer, in list
		liste of  coordinate of G4 for a transcript (start and end)
	    borderTranscript : couple of integer, in list
		liste of couple coordinate of border for a transcript
		borderInferior== border near to 5'
		borderSuperior== border near to 3'
	    Returns
	    -------
	    inTranscript: boolean
		answer if G4 detected is on this transcript
	"""
	startG4=int(coordG4[0])
	endG4=int(coordG4[1])
	borderInferior=int(borderTranscript[0])
	borderSuperior=int(borderTranscript[1])
	inTranscript=False
	#~ print borderInferior, borderSuperior, coordG4
	if (strand == str(1)): ## gene positif
		if (int(startG4) >= int(borderInferior) and int(startG4) <= int(borderSuperior) and int(endG4) >= int(borderInferior) and int(endG4) <= int(borderSuperior)):
			inTranscript=True
	else:## gene negatif
		if (int(startG4) <= int(borderInferior) and int(startG4) >= int(borderSuperior) and int(endG4) <= int(borderInferior) and int(endG4) >= int(borderSuperior)):
			inTranscript=True
	return inTranscript
######################################################################################################################################################
def GetLocalisationPositif(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList):
	""" Give localisation of G4 region in the transcript positif
	    Parameters
	    ----------
	    BiotypeByTranscript : dico
		where key = ID of transcript, value = biotype of this transcript
	    ProteinCoding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    transcriptId : string
		identifiant of the transcript analyse
	    borderTranscript : couple of integer, in list
		liste of couple coordinate of border for a transcript
		borderInferior== border near to 5'
		borderSuperior== border near to 3'
	    coordG4 : couple of integer, in list
		liste of  coordinate of G4 for a transcript (start and end)
	    end5 : string
		end of 5' region, string and not integer because can be empty ('') when transcript doesn't have 5'region
	    start3 :  string
		start of 3' region, string and not integer because can be empty ('') when transcript doesn't have 3'region
	    exonList : list 
 		list of exons contained in this transcript
	    intronList : list
		list of introns contained in this transcript   
	    Returns
	    -------
	    localisation: string
		localisation of G4 region in the transcript positif
	"""
	startG4=int(coordG4[0]) 
	endG4=int(coordG4[1])
	borderInferior=int(borderTranscript[0])
	borderSuperior=int(borderTranscript[1]) 
	localisation='NAP'
	biotypeTranscript= BiotypeByTranscript.get(transcriptId)
	if (biotypeTranscript in ProteinCoding):
		codant='CDS'
	else:

		codant='ExonNC'
	if (end5 != "" and startG4 <= int(end5)):	# if G4 is positioned in the 5' region
		if (endG4 <= int(end5)): # if G4 is positioned totally in the 5' region
			localisation="5"
		else:	# if G4 is positioned in the 5' region with overlapping
			localisation="junction_5_"+codant
	elif (start3 != "" and endG4 >= int(start3)):	# if G4 is positioned in the 3' region
		if (startG4 >= int(start3)):	 # if G4 is positioned totally in the 3' region
			localisation="3"
		else:	# if G4 is positioned in the 3' region with overlapping
			localisation="junction_"+codant+"_3" 
	else:
		for couple_exon in exonList: # if G4 is positioned in the exonic region, for each exon
			if not (couple_exon==""): # if transcript contain exon
				startExon=int(couple_exon.split("-")[0])
				endExon=int(couple_exon.split("-")[1])
				if (startG4 >= startExon and startG4 <= endExon):
					if (endG4 <= endExon):
						localisation=codant
					else:
						if (codant =='CDS'):
							localisation="junction_"+codant+"_Intron"
						else:
							localisation="junction_"+codant+"_IntronNC"
		for couple_Intron in intronList:	# if G4 is positioned in the intronic region, for each intron
			if not (couple_Intron==""):	# if transcript contain intron
				startIntron=int(couple_Intron.split("-")[0])
				endIntron=int(couple_Intron.split("-")[1])
				if (startG4 >= startIntron and startG4 <= endIntron):
					if (endG4 <= endIntron):
						if (codant=='CDS'):
							localisation="Intron"
						else:
							localisation="IntronNC"
					else:
						if (codant=='CDS'):
							localisation="junction_Intron_"+codant
						else:
							localisation="junction_IntronNC_"+codant
	if ((localisation is not None) or (localisation != 'NA')):
		return localisation

######################################################################################################################################################
def GetLocalisationNegatif(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList ):
	""" Give localisation of G4 region in the transcript negatif
	    Parameters
	    ----------
	    BiotypeByTranscript : dico
		where key = ID of transcript, value = biotype of this transcript
	    ProteinCoding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    transcriptId : string
		identifiant of the transcript analyse
	    borderTranscript : couple of integer, in list
		liste of couple coordinate of border for a transcript
		borderInferior== border near to 5'
		borderSuperior== border near to 3'
	    coordG4 : couple of integer, in list
		liste of  coordinate of G4 for a transcript (start and end)
	    end5 : string
		end of 5' region, string and not integer because can be empty ('') when transcript doesn't have 5'region
	    start3 :  string
		start of 3' region, string and not integer because can be empty ('') when transcript doesn't have 3'region
	    exonList : list 
 		list of exons contained in this transcript
	    intronList : list
		list of introns contained in this transcript   
	    Returns
	    -------
	    localisation: string
		localisation of G4 region in the transcript negatif
	"""
	startG4=int(coordG4[0]) #
	endG4=int(coordG4[1])
	borderInferior=int(borderTranscript[0])
	borderSuperior=int(borderTranscript[1]) 
	localisation='NAN'
	biotypeTranscript= BiotypeByTranscript.get(transcriptId)
	if (biotypeTranscript in ProteinCoding):
		codant='CDS'
	else:

		codant='ExonNC'
	if (end5 != "" and endG4 >= int(end5)):	# if G4 is positioned in the 5' region
		if (startG4 >= int(end5)): # if G4 is positioned totally in the 5' region
			localisation="5"
		else:	# if G4 is positioned in the 5' region with overlapping
			localisation="junction_5_"+codant
	elif (start3 != "" and startG4 <= int(start3)):	# if G4 is positioned in the 3' region
		if (endG4 <= int(start3)):	 # if G4 is positioned totally in the 3' region
			localisation="3"
		else:	# if G4 is positioned in the 3' region with overlapping
			localisation="junction_"+codant+"_3" 
	else:
		for couple_exon in exonList: # if G4 is positioned in the exonic region, for each exon
			if not (couple_exon==""): # if transcript contain exon
				startExon=int(couple_exon.split("-")[0])
				endExon=int(couple_exon.split("-")[1])
				if (endG4 <= startExon and endG4 >= endExon):
					if (startG4 >= endExon):
						localisation=codant
					else:
						if (codant =='CDS'):
							localisation="junction_"+codant+"_Intron"
						else:
							localisation="junction_"+codant+"_IntronNC"
		for couple_Intron in intronList:	# if G4 is positioned in the intronic region, for each intron
			if not (couple_Intron==""):	# if transcript contain intron
				startIntron=int(couple_Intron.split("-")[0])
				endIntron=int(couple_Intron.split("-")[1])
				if (endG4 <= startIntron and endG4 >= endIntron):
					if (startG4 >= endIntron):
						if (codant=='CDS'):
							localisation="Intron"
						else:
							localisation="IntronNC"
					else:
						if (codant=='CDS'):
							localisation="junction_Intron_"+codant
						else:
							localisation="junction_IntronNC_"+codant
	if ((localisation is not None) or (localisation != 'NA')):
		return localisation
######################################################################################################################################################
def GetAnnotationTranscript(filename,ProteinCoding,BiotypeByTranscript):
	""" This fonction is to define which transcript contain a good annotation or not from Ensembl
	    Creation of a dico, where key = id of the transcript and value= bolean value (False if the annotation isn't good)
	    Parameters
	    ----------
	    filename : string
		name of the file which contain informatiob fr each transcript (except the transcript's biotype)
	    ProteinCoding : ensemble of string, in list
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
		biotypeTranscript=BiotypeByTranscript.get(transcriptId)	## get transcriptBiotype of this TRanscriptId
		if (biotypeTranscript not in ProteinCoding): # if transcript not a protein coding
			if (start5!='' or end5!='' or start3!='' or end3!='' ): # but if transcript has a 5'UTR or an 3' UTR
				answer=False # has not a good annotation in Ensembl
		AnnotationTranscript[transcriptId]=answer
	return AnnotationTranscript
######################################################################################################################################################
def GetLocalisationG4InJunction(BiotypeByTranscript,ProteinCoding,transcriptId):
	localisation='NAJ'
	biotypeTranscript= BiotypeByTranscript.get(transcriptId)
	if (biotypeTranscript in ProteinCoding):
		codant='CDS'
	else:

		codant='ExonNC'
	return "junction_"+codant+"_"+codant
	
######################################################################################################################################################
def GetLocalisationG4InTranscript(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList, strand):
	""" This fonction is to define which fonction will be used for determinate the localisation in transcript. 
		If the transcript has a strand positif, it's will be the fonction "GetLocalisationPositif".
		If the transcript has a strand negatif, it's will be the fonction "GetLocalisationNegatif".
	    ----------
	    BiotypeByTranscript : dico
		dictionary, where key =id of transcript and value = biotype of this transcript
	    ProteinCoding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    transcriptId : string
		identifiant of the transcript analyse
	    borderTranscript : couple of integer, in list
				liste of couple coordinate of border for a transcript
				borderInferior== border near to 5'
				borderSuperior== border near to 3'
	    end5 : string
		end of 5' region, string and not integer because can be empty ('') when transcript doesn't have 5'region
	    start3 :  string
		start of 3' region, string and not integer because can be empty ('') when transcript doesn't have 3'region
	    exonList : list 
 		list of exons contained in this transcript
	    intronList : list
		list of introns contained in this transcript
	    strand : string
		representation of the gene sens (1 for gene positif, -1 for gene negatif)
	    Returns
	    -------
	    localisation: string
		localisation of G4 region in the transcript
	"""
	if (strand == str(1)): # if gene positif
		localisation=GetLocalisationPositif(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList)
	else: # if gene negatif
		localisation=GetLocalisationNegatif(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList)
	return localisation
######################################################################################################################################################	
def AddG4InTranscriptome(G4InTranscript,transcriptId, descriptionG4,informationsOfG4,localisationInTranscript, biotypeTranscript):
	""" This fonction add for each G4 founded the informations issu from G4 screener (cGcC, g4H, sequence and G4NN) to the localisation of G4 and the biotype of transcript
	    ----------
	    G4InTranscript : dico
		dictionary, where key =id of transcript+ description of G4 and value = all info of this G4
	    transcriptId : string
		identifiant of the transcript analyse
	    descriptionG4 : string
		ensemble of information for describe the G4 by it position (chromosome:startG4-endG4|strand)
	    informationsOfG4 : list
		list of informations obtained by G4-screener (cGcC, g4H, sequence and G4NN )
	    localisationInTranscript :  string
		localisation of G4 on this transcript
	    biotypeTranscript : string
 		biotype of this this transcript
	    Returns
	    -------
	    G4InTranscript : dico
		dictionary, where key =id of transcript+ description of G4 and value = all info of this G4
	"""
	value = [] # list empty of value for a g4 by gene
	for information in informationsOfG4: # add information already obtained by G4 screener (cGcC, g4H, sequence and G4NN)
		value.append(information)
	value.append(localisationInTranscript) # add localisation of this G4
	value.append(biotypeTranscript) # add biotype of this G4 for this transcript
	headerG4=transcriptId+'|'+descriptionG4 # header uniq by G4 in transcript
	if (G4InTranscript.has_key(headerG4) == False): # if the transcript doesn't be traited before
		G4InTranscript[headerG4]=value # add this G4 in the dictionary
	return G4InTranscript	# return the dictionary
######################################################################################################################################################
def AddG4InGenome(G4InGenome, geneId, descriptionG4, localisationInTranscript):
	""" This fonction add for each G4 founded in the gene all localisation possible
	    ----------
	    G4InGenome : dico
		dictionary, where key =id of gene+ description of G4 and value = all localisation of this G4
	    geneId : string
		identifiant of the gene analyse
	    descriptionG4 : string
		ensemble of information for describe the G4 by it position (chromosome:startG4-endG4|strand)
	    informationsOfG4 : list
		list of informations obtained by G4-screener (cGcC, g4H, sequence and G4NN )
	    localisationInTranscript :  string
		localisation of G4 on this transcript
	    Returns
	    -------
	    G4InTranscript : dico
		dictionary, where key =id of transcript+ description of G4 and value = all info of this G4
	"""
	headerG4=geneId+':'+descriptionG4  # header uniq by G4 in gene
	if (G4InGenome.has_key(headerG4) == False ): # if the gene doesn't be traited before
		G4InGenome[headerG4]=localisationInTranscript 	# traited it
	else:	# if the gene already traited before
		localisationInGene=G4InGenome.get(headerG4) # get the localisation already known for this G4
		if (localisationInGene is not None and localisationInTranscript not in localisationInGene):
			G4InGenome[headerG4]=localisationInGene+';'+localisationInTranscript # add this localisation
	return G4InGenome
######################################################################################################################################################
def AddTranscriptPerG4(TranscriptPerG4, descriptionG4,transcriptId): 

	if (TranscriptPerG4.has_key(descriptionG4) == False): 
		TranscriptPerG4[descriptionG4] = transcriptId
	else:
		if bool(TranscriptPerG4): ## dico already with info
			lastTranscript= TranscriptPerG4.get(descriptionG4)
			TranscriptPerG4[descriptionG4] = lastTranscript +"-"+transcriptId
	return TranscriptPerG4
######################################################################################################################################################
def JuncionIsInTranscript(coordG4,intronList):
	""" Give boolean, if a junction where was found a G4 is located in a transcript
	    Parameters
	    ----------
	    coordG4 : couple of integer, in list
		liste of  coordinate of G4 for a transcript (start and end)
	     intronList : list
		list of introns contained in this transcript
	    Returns
	    -------
	    inTranscript: boolean
		answer if G4 detected on junction is on this transcript
	"""
	startG4=int(coordG4[0]) 
	endG4=int(coordG4[1])
	inTranscript=False
	for couple_Intron in intronList: # for each couple of intron of this transcript
		if not (couple_Intron==""):
			endIntron=int(couple_Intron.split("-")[1])	# end of this intron
			if (startG4 <= endIntron and endIntron<= endG4):
				inTranscript=True
	return inTranscript
######################################################################################################################################################
def GetlisteG4InGene(G4Detected, listeG4InGene):
	""" This fonction is to create a dictionary which represent for each gene (geneId) the list of G4 detected for it, depend to the dictionary 
	use in entry (G4 detected in sequence entire of the gene or in the section junction created)
	    ----------
	    G4Detected : dico
		dictionary, where key =description of G4 (geneId|startG4|endG4|strand) and value = biotype of this transcript
	    listeG4InGene : dico
		dictionary, where key =id of gene and value = list of G4 (description of G4) found in this gene
	    Returns
	    -------
	    listeG4InGene: dico
		dictionary, where key =id of gene and value = list of G4 (description of G4) found in this gene
	"""
	for descriptionG4, informationsOfG4 in G4Detected.items():
		geneID=descriptionG4.split('|')[0]
		if (listeG4InGene.has_key(geneID) == False):	
			listeG4=[]
		else:
			listeG4=listeG4InGene.get(geneID)
		listeG4.append(descriptionG4)
		listeG4InGene[geneID]=listeG4
	return listeG4InGene

######################################################################################################################################################
def ExtractionG4InTranscript(directory, specie, chromosome, G4InTranscript):
	output= open(directory+"/"+specie+"_chr"+chromosome+"_G4InTranscript.txt","w") ## file opening
	output.write("InfoG4ByTranscript\tcGcC\tG4Hunter\tsequenceG4\tG4NN\tlocalisation\ttranscriptBiotype\n")
	for key,value in G4InTranscript.items():
		if None not in value: # because some transcriptID from ensembl donMt contain info of biotype (as ENST00000604369)
			output.write(key+"\t"+'\t'.join(value)+"\n")

def ExtractionG4InGenome(directory, specie, chromosome, G4InGenome):
	output= open(directory+"/"+specie+"_chr"+chromosome+"_G4InGenome.txt","w") ## file opening
	output.write("InfoG4ByGene\tLocalisation(s)\n")
	for key,value in G4InGenome.items():
		output.write(key+"\t"+value+"\n")

def ExtractionTranscriptPerG4(directory, specie, chromosome, TranscriptPerG4):
	output= open(directory+"/"+specie+"_chr"+chromosome+"_TranscriptPerG4.txt","w") ## file opening
	output.write("CoordonnéesG4\tTranscript(s)\n")
	for key,value in TranscriptPerG4.items():
		output.write(key+"\t"+value+"\n")
######################################################################################################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	parser.add_argument ('-pG4', '--pathG4', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/mouseEsssai/')
	parser.add_argument ('-pBtp', '--pathBiotype', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/mouse/transcriptType')
	parser.add_argument ('-o', '--outPut', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/mouse/G4Screener')
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
	pathG4=arg.pathG4	# directory which contain all the directory chromosome
	pathBiotype=arg.pathBiotype # directory which contain all the biotype
	outPut=arg.outPut # directory used for output
	CHROMOSOME=arg.CHROMOSOME	# chromosome to analyze
	specie=arg.specie	# specie to analyse
	THRESHOLD_G4H=arg.THRESHOLD_G4H	# threshold use to discriminate the score G4H (litterature = 0.9)
	THRESHOLD_CGCC=arg.THRESHOLD_CGCC	# threshold use to discriminate the score G4H (litterature = 4.5)
	THRESHOLD_G4NN=arg.THRESHOLD_G4NN	# threshold use to discriminate the score G4NN (litterature = 0.5)
	EXTENSION=arg.EXTENSION	# EXTENSION use for create the junction exon_exon previously
	WINDOW=arg.WINDOW	# size of window use in G4 screener
	STEP=arg.STEP	# size of step use in G4 screener



	ProteinCoding=['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene', 'IG_M_gene', 'IG_V_gene', 'IG_Z_gene', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay', 'protein_coding', 'TR_C_gene', 'TR_D_gene', 'TR_gene', 'TR_J_gene', 'TR_V_gene']

	BiotypeByTranscript={}	# dictionary of strand for each transcript
	G4DetectedInGene={}
	G4DetectedInJunction={}
	StrandByGene={}
	listeG4InGeneEntire={}
	listeG4InGeneJunction={}

	directory=pathG4+'/chr'+CHROMOSOME+'_csv'	# variable directory which contain the data for this chromosome
	index=directory+'/'+specie+'_transcript_unspliced_chr'+CHROMOSOME+'_Index.txt'	# file which contain info by transcript for this chromosome
	indexBiotypeTranscript=pathBiotype+'/transcriptType_chr'+CHROMOSOME+'.txt'	# file which contain biotype of transcript for this chromosome

	BiotypeByTranscript=CreateDictionaryBiotypeByTranscript (indexBiotypeTranscript) # dictionary of biotype for each transcript
	StrandByGene=CreateDictionaryStrandByGene (index)	# dictionary of strand for each transcript
	
	AnnotationTranscript={}
	AnnotationTranscript=GetAnnotationTranscript(index,ProteinCoding,BiotypeByTranscript) # annotation transcript
	
	

	
	
	
	#####Creation dictionary G4_detected_gene and G4_detected_junction , for G4 found in gene and junction with scores > thresolds
	for path, dirs, files in os.walk(directory): # for each element of the directory to passed
		for filename in files: # for each files
			inputfile=directory+'/'+filename
			parametersTool=[THRESHOLD_CGCC,THRESHOLD_G4H,THRESHOLD_G4NN, WINDOW, STEP]
			##EXTRACTION OF G4 BY SCORE G4NN > THRESHOLD IN A DICO
			if ('gene_unspliced' in filename and '.csv' in filename ): ## for G4 in gene	
				G4DetectedInGene=ReturnG4InGene(G4DetectedInGene, inputfile, parametersTool, StrandByGene)
			elif ('transcript_unspliced' in filename and '.csv' in filename): ## for G4 in junction CDS-CDS --> from splicing
				G4DetectedInJunction=ReturnG4InJunction(G4DetectedInJunction, inputfile, parametersTool, EXTENSION, StrandByGene)


	
	listeG4InGeneEntire=GetlisteG4InGene(G4DetectedInGene, listeG4InGeneEntire)
	listeG4InGeneJunction=GetlisteG4InGene(G4DetectedInJunction, listeG4InGeneJunction)
	


	G4InTranscript={}
	G4InGenome={}
	TranscriptPerG4={}
	
	
	############################################################################## g4 entire
	inputfile= open(index,"r") # file opening for reading
	for line in inputfile: # for each transcript
			words=line.split('|')
			transcriptId=words[0].rstrip() 
			geneId=words[1].rstrip() 
			chromosome=words[2].rstrip() 
			strand=words[3].rstrip() 
			biotypeGene=words[4].rstrip() 
			exonList=words[5].rstrip().split(";") ## transform in array
			intronList=words[6].rstrip().split(";") ## transform in array
			### resolution probleme start5 and end5 for gene negatif
			if (strand == str(1)):
				start5=words[7].rstrip()
				end5=words[8].rstrip()
				start3=words[9].rstrip()
				end3=words[10].rstrip()
			else:

				end5=words[7].rstrip()
				start5=words[8].rstrip() 
				end3=words[9].rstrip()
				start3=words[10].rstrip()
			
			biotypeTranscript=BiotypeByTranscript.get(transcriptId) # get biotype of this transcript
			borderTranscript=BorderOfTranscript(start5 , end5 , start3 , end3, exonList, intronList, strand)
			if (listeG4InGeneEntire.has_key(geneId) == False): # if gene not contain G4 detected
				continue
			else: # if gene contain G4 detected
				listeG4InGene=listeG4InGeneEntire.get(geneId)
				for G4InGene in listeG4InGene:
					startG4=int(G4InGene.split('|')[1])
					endG4=int(G4InGene.split('|')[2])
					coordG4=[startG4,endG4]
					g4inTranscript=G4IsInTranscript(strand, coordG4, borderTranscript)
					annotationTranscript=AnnotationTranscript.get(transcriptId)
					
					if (g4inTranscript == True and annotationTranscript ==True):
						informationsOfG4=list(G4DetectedInGene.get(G4InGene))
						localisationInTranscript=GetLocalisationG4InTranscript(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList, strand)	
						descriptionG4=chromosome+':'+str(startG4)+'-'+str(endG4)+'|'+strand
						G4InTranscript=AddG4InTranscriptome(G4InTranscript,transcriptId, descriptionG4,informationsOfG4,localisationInTranscript, biotypeTranscript)
						G4InGenome=AddG4InGenome(G4InGenome, geneId, descriptionG4, localisationInTranscript)
						TranscriptPerG4=AddTranscriptPerG4(TranscriptPerG4, descriptionG4,transcriptId)
					#	if (localisationInTranscript=='NAP' or localisationInTranscript=='NAN' ):
						#~ print localisationInTranscript
							#~ print  localisationInTranscript, descriptionG4, exonList, start5, end5, start3, end3, strand
							#~ print '-------------------------------'

	############################################################################## g4 junction
	inputfile= open(index,"r") # file opening for reading
	for line in inputfile: # for each transcript
			words=line.split('|')
			transcriptId=words[0].rstrip() 
			geneId=words[1].rstrip() 
			chromosome=words[2].rstrip() 
			strand=words[3].rstrip() 
			biotypeGene=words[4].rstrip() 
			exonList=words[5].rstrip().split(";") ## transform in array
			intronList=words[6].rstrip().split(";") ## transform in array
			### resolution probleme start5 and end5 for gene negatif
			if (strand == str(1)):
				start5=words[7].rstrip()
				end5=words[8].rstrip()
				start3=words[9].rstrip()
				end3=words[10].rstrip()
			else:

				end5=words[7].rstrip()
				start5=words[8].rstrip() 
				end3=words[9].rstrip()
				start3=words[10].rstrip()
			biotypeTranscript=BiotypeByTranscript.get(transcriptId) # get biotype of this transcript
			if (listeG4InGeneJunction.has_key(geneId) == False): # if gene not contain G4 detected
				continue
			else:	# if gene contain G4 detected
				listeG4InGene=listeG4InGeneJunction.get(geneId)
				for G4InGene in listeG4InGene:
					startG4=int(G4InGene.split("|")[1])
					endG4=int(G4InGene.split("|")[2])
					coordG4=[startG4,endG4]
					junctionInTranscript=JuncionIsInTranscript(coordG4,intronList)
					annotationTranscript=AnnotationTranscript.get(transcriptId)
					if (junctionInTranscript == True and annotationTranscript ==True):
						informationsOfG4=list(G4DetectedInJunction.get(G4InGene))
						localisationInTranscript=GetLocalisationG4InJunction(BiotypeByTranscript,ProteinCoding,transcriptId)
					
						descriptionG4=chromosome+':'+str(startG4)+'-'+str(endG4)+'|'+strand
						G4InTranscript=AddG4InTranscriptome(G4InTranscript,transcriptId, descriptionG4,informationsOfG4,localisationInTranscript, biotypeTranscript)
						G4InGenome=AddG4InGenome(G4InGenome, geneId, descriptionG4, localisationInTranscript)
						TranscriptPerG4=AddTranscriptPerG4(TranscriptPerG4, descriptionG4,transcriptId)
				
	
	
	ExtractionG4InTranscript(outPut, specie, CHROMOSOME, G4InTranscript)	
	ExtractionG4InGenome(outPut, specie, CHROMOSOME, G4InGenome)	
	ExtractionTranscriptPerG4(outPut, specie, CHROMOSOME, TranscriptPerG4)
	

main()
	
