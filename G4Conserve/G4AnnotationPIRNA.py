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

import sys
import os
import re
import argparse


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
	length = 0
	if (positionA and positionB) != "":
		# some transcript doesn't have 5' 
		length = ( abs( int(positionA) - int(positionB) ) ) + 1
	return length

def ReturnG4InGene(G4DetectedInGene, inputfile, parametersTool):
	""" Add informations in dictionary with informations of region G4 detected (score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		regions of G4 discovered in the genes

	    Parameters
	    ----------
	    G4DetectedInGene : dictionary 
		dictionary will contain the information of G4 (score cgCc, scote G4Hunter, sequence and score G4NN) for all the region G4 detected
	    inputfile : string
		name of file which contain the liste of informations for each piRNA 
	    parametersTool : list of value use by G4 screener
		parametersTool[0]=float constante,threshold use to discriminate the score cGcC
		parametersTool[1]= float constante,threshold use to discriminate the score G4H
		parametersTool[2]= float constante,threshold use to discriminate the score G4NN
		parametersTool[3]= integer constante,size of window use in G4 screener	 
		parametersTool[4]=integer constante,size of step use in G4 screener
	    Returns
	    -------
	    G4DetectedInGene
		dictionary with informations of region G4 detected (score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		regions of G4 discovered in the genes
	"""
	rnaType='piRNA'
	localisation='ExonNC'
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
			infopiRNA=words[1].rstrip()
			cGcC=float(words[2].rstrip()) #  score cGcC for this window
			g4H=float(words[3].rstrip()) #  score G4Hunter for this window 
			sequence=words[4].rstrip()	# sequence for this window
			startWindow=int(words[5].rstrip())	# start window from G4screener
			endWindow=int(words[6].rstrip()) ## end window from G4screener
			g4NN=float(words[7].rstrip()) # score G4NN for this window
			if (cGcC>=THRESHOLD_CGCC and g4H>=THRESHOLD_G4H and g4NN>=THRESHOLD_G4NN): # if window contain value > several threshold
				G4DetectedInGene[infopiRNA]=str(cGcC), str(g4H),sequence , str(g4NN), localisation, rnaType
			else: # if it's not the first windows above the thresholds
				continue

	return G4DetectedInGene
######################################################################################################################################################
def getLengthpiRNA (inputfile):
	""" Count length of all piRNA present in the file input
	    Parameters
	    ----------
	    inputfile : string
		name of file which contain the liste of informations for each piRNA 
	    Returns
	    -------
	    length
		int of length of all piRNA present in the file 
	"""
	length=0
	inputfile= open(inputfile,"r")	# file opening for reading
	for line in inputfile:	# for each file csv create by G4 screener
		if (re.search('^[0-9]', line)): # if the line is not the header of the file 
			words=line.split('\t') # parsing by the separator, here '\t' (file csv)
			numRow=words[0].rstrip() #numero of row for this for this window
			infopiRNA=words[1].rstrip()
			cGcC=float(words[2].rstrip()) #  score cGcC for this window
			g4H=float(words[3].rstrip()) #  score G4Hunter for this window 
			sequence=words[4].rstrip()	# sequence for this window
			startWindow=int(words[5].rstrip())	# start window from G4screener
			endWindow=int(words[6].rstrip()) ## end window from G4screener
			g4NN=float(words[7].rstrip()) # score G4NN for this window
			length=length+len(sequence)
	return length
######################################################################################################################################################
def getNbrpiRNAWithPG4(dictionary): 
	""" Count number piRNA which contain the PG4s
	    Parameters
	    ----------
	    dictionary : 
		dictionary with informations of region G4 detected (score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		regions of G4 discovered in the genes
	    Returns
	    -------
	    nbrpiRNAWithPG4
		int of number piRNA which contain the PG4s 
	"""
	return list(dictionary) # list of PG4 detected in the piRNA (key)
######################################################################################################################################################
def getNbrpiRNA(inputfile): 
	""" Count length of all piRNA present in the file input
	    Parameters
	    ----------
	    inputfile : string
		name of file which contain the liste of informations for each piRNA 
	    Returns
	    -------
	   nbrpiRNAWithPG4
		int of number piRNA test 
	"""
	listepiRNA=[]
	inputfile= open(inputfile,"r")	# file opening for reading
	for line in inputfile:	# for each file csv create by G4 screener
		if (re.search('^[0-9]', line)): # if the line is not the header of the file 
			words=line.split('\t') # parsing by the separator, here '\t' (file csv)
			numRow=words[0].rstrip() #numero of row for this for this window
			infopiRNA=words[1].rstrip()
			cGcC=float(words[2].rstrip()) #  score cGcC for this window
			g4H=float(words[3].rstrip()) #  score G4Hunter for this window 
			sequence=words[4].rstrip()	# sequence for this window
			startWindow=int(words[5].rstrip())	# start window from G4screener
			endWindow=int(words[6].rstrip()) ## end window from G4screener
			g4NN=float(words[7].rstrip()) # score G4NN for this window
			listepiRNA.append(infopiRNA)
	return listepiRNA
		
		

######################################################################################################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	parser.add_argument ('-i', '--inputfile', default = '/home/local/USHERBROOKE/bels2814/Documents/These/piRNA/piR_human_v1.0.csv')
	parser.add_argument ('-specie', '--specie', default = 'HS')
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
	inputfile=arg.inputfile	# directory which contain all the directory chromosome
	THRESHOLD_G4H=arg.THRESHOLD_G4H	# threshold use to discriminate the score G4H (litterature = 0.9)
	THRESHOLD_CGCC=arg.THRESHOLD_CGCC	# threshold use to discriminate the score G4H (litterature = 4.5)
	THRESHOLD_G4NN=arg.THRESHOLD_G4NN	# threshold use to discriminate the score G4NN (litterature = 0.5)
	EXTENSION=arg.EXTENSION	# EXTENSION use for create the junction exon_exon previously
	WINDOW=arg.WINDOW	# size of window use in G4 screener
	STEP=arg.STEP	# size of step use in G4 screener


	G4DetectedInGene={}
	length=0
	piRNAWithPG4=[]
	piRNA=[]


	parametersTool=[THRESHOLD_CGCC,THRESHOLD_G4H,THRESHOLD_G4NN, WINDOW, STEP]
		
	

	G4DetectedInGene=ReturnG4InGene(G4DetectedInGene, inputfile, parametersTool) ### detection of PG4 in piRNA case

	length=getLengthpiRNA (inputfile) # length of all piRNA
	piRNAWithPG4=getNbrpiRNAWithPG4 (G4DetectedInGene) ## list of piRNA with at least one PG4
	piRNA=getNbrpiRNA(inputfile) ## list of piRNA tested
		
	
	print '\nInfo PG4 detected in piRNA:'
	print '\tName File  :'+inputfile
	if ( not G4DetectedInGene):
		print '		NO PG4 DETECTED'

	print '\nNumber PG4 detected in piRNA:	'+str(len(G4DetectedInGene))
	print 'Number piRNA with PG4:	'+str(len(piRNAWithPG4))
	print 'Number piRNA tested:	'+str(len(piRNA))
	print 'Length of all piRNA (nt):	'+str(length)
	print 'Densité in Kb :	'+str(round(len(G4DetectedInGene)/ float(length) * 1000,2))
	print 'Percent piRNA with Kb:	'+str(round(len(G4DetectedInGene)/ float(len(piRNA)) * 100,2))
	
main()
	
