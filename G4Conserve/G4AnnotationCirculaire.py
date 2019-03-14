#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
Copyright Université de Sherbrooke, Département de Biochimie et Département d'Informatique

contact: sarah.belhamiti@usherbrooke.ca

This software is a computer program whose annote the G4 region of one chromosome of one specie.


---------------------------------------------------

``G4AnnotationCirculaire`` **module description**:

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
def ReturnG4InCircu(G4Detected, inputfile, parametersTool): ## for G4 in gene
	""" Add informations in dictionary with informations of region G4 detected (score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		regions of G4 discovered in the genes

	    Parameters
	    ----------
	    G4Detected : dictionary 
		dictionary will contain the information of G4 (score cgCc, scote G4Hunter, sequence and score G4NN) for all the region G4 detected
	    inputfile : string
		name of file which contain the liste of informations for each circular 
	    parametersTool : list of value use by G4 screener
		parametersTool[0]=float constante,threshold use to discriminate the score cGcC
		parametersTool[1]= float constante,threshold use to discriminate the score G4H
		parametersTool[2]= float constante,threshold use to discriminate the score G4NN
		parametersTool[3]= integer constante,size of window use in G4 screener	 
		parametersTool[4]=integer constante,size of step use in G4 screener
	    Returns
	    -------
	    G4Detected
		dictionary with informations of region G4 detected (score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		regions of G4 discovered in the genes
	"""
	'''
	rnaType='circu'
	localisation='ExonNC'
	THRESHOLD_CGCC=parametersTool[0]
	THRESHOLD_G4H=parametersTool[1]
	THRESHOLD_G4NN=parametersTool[2]
	WINDOW=parametersTool[3]
	STEP=parametersTool[4]
	oldPassed=False
	passed= False
	overtakingSequence=False
	first=False
	tmp=[]
	amont=''
	aval=''
	descriptionOverThreshold=''
	inputfile= open(inputfile,"r")	# file opening for reading
	for line in inputfile:	# for each file csv create by G4 screener
		if (re.search('^[0-9]', line)): # if the line is not the header of the file 
			words=line.split('\t') # parsing by the separator, here '\t' (file csv)
			if (len(words) == 8):
				numRow=words[0].rstrip() #numero of row for this for this window
				infoCircu=words[1].rstrip()
				cGcC=float(words[2].rstrip()) #  score cGcC for this window
				g4H=float(words[3].rstrip()) #  score G4Hunter for this window 
				sequence=words[4].rstrip()	# sequence for this window
				startWindow=int(words[5].rstrip())	# start window from G4screener
				endWindow=int(words[6].rstrip()) ## end window from G4screener
				g4NN=float(words[7].rstrip()) # score G4NN for this window
				idCircu=infoCircu.split('|')[0] # id sequence circular
				chromosome=infoCircu.split('|')[1].split(':')[0] # chromosome
				startBorder=infoCircu.split('|')[1].split(':')[1].split('-')[0]
				endBorder=infoCircu.split('|')[1].split(':')[1].split('-')[1].replace("+","")
				strand=infoCircu.split('|')[1].split(':')[1][-1]
				description=chromosome+'|'+str(startBorder)+'-'+str(endBorder)

				if (tmp and descriptionOverThreshold !=description):
					print '333333333333', headerG4
					if (strand == '+'):
						headerG4=str(tmp[0])+"|"+str(tmp[1])+':'+str(tmp[2])+"-"+str(tmp[3])+"|"+'1'
					else:
						headerG4=str(tmp[0])+"|"+str(tmp[1])+':'+str(tmp[2])+"-"+str(tmp[3])+"|"+'-1'
					localisation='ExonNC'
					if (G4Detected.has_key(headerG4) == False and strand != None): 
						G4Detected[headerG4]=str(tmp[4]), str(tmp[5]),tmp[6] , str(tmp[7]), localisation, rnaType
					oldPassed=False
					passed= False
					overtakingSequence=False
					first=False
					tmp=[]
					amont=''
					aval=''			
				
				if (cGcC>=THRESHOLD_CGCC and g4H>=THRESHOLD_G4H and g4NN>=THRESHOLD_G4NN and overtakingSequence == False): # if window contain value > several threshold
					passed=True # passage over the differents threshold for this window
					if (int(startWindow) == 1): # if window is the first of the sequence
						first=True
					else: # if window is not the first of the sequence
						first=False
					if (oldPassed != passed): # if it's the first windows (beginning if the passage
						listeCGcC=[] # list which will contain the values of score cgCG for windows over the THRESHOLD_CGCC
						listeG4H=[]  # list which will contain the values of score G4H for windows over the THRESHOLD_G4H
						listeG4NN=[]  # list which will contain the values of score G4NN for windows over the THRESHOLD_G4NN
						descriptionOverThreshold=description # assignation of description for this variable
						sequenceG4=sequence
						oldPassed=passed # assignation, we pass over the thresholds
						if (strand == '+'): # if gene positif
							startG4=startWindow # startG4 same as startWindow of G4 screener
							endG4=endWindow # endG4 same as endWindow of G4 screener
						else:	# if gene negatif
							startG4=int(endBorder)-int(startWindow)+1 # calcul of startG4
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
					
						if (strand == '+'): # if gene positif
							endG4=endWindow # endG4 moves from the end of the first window,  same as endWindow of G4 screener
						else: # if gene negatif
							endG4=int(startG4)-(len(sequenceG4))+1 # endG4 moves from the end of the first window, # calcul of endG4
						lastRow=numRow # numRow of the last windows study over the threshols
						listeCGcC.append(cGcC) # add score cGcC for the windows over the threshold
						listeG4Hunter.append(g4H) # add score G4H for the windows over the threshold
						listeG4NN.append(g4NN) # add score G4NN for the windows over the threshold
					#### test condition si depassement de la sequence
					if (overtakingSequence == False): # if border sequence not depass
						if (strand == '+'): # if gene positif
							if (int(startG4) > int(startBorder) and int(startG4) < int(endBorder) and int(endG4) > int(endBorder)):
								overtakingSequence=True # depassement
								
					
						else:
							if (int(startG4) > int(startBorder) and int(startG4) < int(endBorder) and int(endG4) < int(startBorder)):
								overtakingSequence=True # depassement
						
				if (cGcC<THRESHOLD_CGCC or g4H<THRESHOLD_G4H or g4NN<THRESHOLD_G4NN or descriptionOverThreshold != description or test==True):
				# if one of the score is under his threshold or if this windows contain info of an other gene
					passed=False # passed passed became false (pass under the thresolds)
					if (oldPassed != passed ): # last windows before under the thresolds
						meanCGcC=mean(listeCGcC) # mean of the score cGcC will be the score for this region containing G4
						meanG4Hunter=mean(listeG4Hunter) # mean of the score G4Hunter will be the score for this region containing G4
						meanG4NN=mean(listeG4NN) # mean of the score G4NN will be the score for this region containing G4
						oldPassed=passed # assignation, we pass under the thresholds
						if (overtakingSequence==True and test==True): # if depassement 
							if (strand == '+' ):
									endG4=tmp[2]
									amont=sequenceG4[:GetLengthFraction(startBorder,endG4)-1]
									aval=tmp[6][:GetLengthFraction(startBorder,startG4)-1]
									sequenceG4=amont+aval
									meanCGcC=mean([meanCGcC,tmp[4]])
									meanG4Hunter=mean([meanG4Hunter,tmp[5]])
									meanG4NN=mean([meanG4NN,tmp[7]])
									localisation='junction'
									headerG4=idCircu+"|"+chromosome+':'+str(startG4)+"-"+str(endG4)+"|"+'1'
							elif (strand == '-'):
									endG4=tmp[2]
									amont=sequenceG4[:GetLengthFraction(startBorder,startG4)-1]
									aval=tmp[6][:GetLengthFraction(endBorder,startG4)-1]
									sequenceG4=amont+aval	
									meanCGcC=mean([meanCGcC,tmp[4]])
									meanG4Hunter=mean([meanG4Hunter,tmp[5]])
									meanG4NN=mean([meanG4NN,tmp[7]])
									localisation='junction'
									headerG4=idCircu+"|"+chromosome+':'+str(endG4)+"-"+str(startG4)+"|"+'-1'
							ptin
							print '999999999999999', headerG4
							G4Detected[headerG4]=str(meanCGcC), str(meanG4Hunter),sequenceG4 , str(meanG4NN), localisation, rnaType
							oldPassed=False
							passed= False
							overtakingSequence=False
							first=False
							tmp=[]
							amont=''
							aval=''
						else: # not depassing	
							if (first != True ): # if isnt firts windows 
								if (strand == '+'):
									headerG4=idCircu+"|"+chromosome+':'+str(startG4)+"-"+str(endG4)+"|"+'1'
								elif (strand == '-'):
									headerG4=idCircu+"|"+chromosome+':'+str(endG4)+"-"+str(startG4)+"|"+'-1'
								if (G4Detected.has_key(headerG4) == False and strand != None): 
									print '55555555', headerG4
									G4Detected[headerG4]=str(meanCGcC), str(meanG4Hunter),sequenceG4 , str(meanG4NN), localisation, rnaType
								if (descriptionOverThreshold != description):
									oldPassed=False
									passed= False
									overtakingSequence=False
									first=False
									tmp=[]
									amont=''
									aval=''
							else:
		
									tmp=[idCircu,chromosome,endG4,startG4,meanCGcC,meanG4Hunter,sequenceG4,meanG4NN]
									continue

	if (tmp):
		if (strand == '+'):
			headerG4=str(tmp[0])+"|"+str(tmp[1])+':'+str(tmp[2])+"-"+str(tmp[3])+"|"+'1'
		else:
			headerG4=str(tmp[0])+"|"+str(tmp[1])+':'+str(tmp[2])+"-"+str(tmp[3])+"|"+'-1'
		localisation='ExonNC'
		print '2222', headerG4
		if (G4Detected.has_key(headerG4) == False and strand != None): 
			
			G4Detected[headerG4]=str(tmp[4]), str(tmp[5]),tmp[6] , str(tmp[7]), localisation, rnaType
	'''

	rnaType='circu'
	localisation='ExonNC'
	THRESHOLD_CGCC=parametersTool[0]
	THRESHOLD_G4H=parametersTool[1]
	THRESHOLD_G4NN=parametersTool[2]
	WINDOW=parametersTool[3]
	STEP=parametersTool[4]
	oldPassed=False
	passed= False
	overtakingSequence=False
	first=False
	tmp=[]
	amont=''
	aval=''
	condition=False
	descriptionOverThreshold=''
	inputfile= open(inputfile,"r")	# file opening for reading
	for line in inputfile:	# for each file csv create by G4 screener
		if (re.search('^[0-9]', line)): # if the line is not the header of the file 
			words=line.split('\t') # parsing by the separator, here '\t' (file csv)
			if (len(words) == 8):
				numRow=words[0].rstrip() #numero of row for this for this window
				infoCircu=words[1].rstrip()
				cGcC=float(words[2].rstrip()) #  score cGcC for this window
				g4H=float(words[3].rstrip()) #  score G4Hunter for this window 
				sequence=words[4].rstrip()	# sequence for this window
				startWindow=int(words[5].rstrip())	# start window from G4screener
				endWindow=int(words[6].rstrip()) ## end window from G4screener
				g4NN=float(words[7].rstrip()) # score G4NN for this window
				idCircu=infoCircu.split('|')[0] # id sequence circular
				chromosome=infoCircu.split('|')[1].split(':')[0] # chromosome
				startBorder=infoCircu.split('|')[1].split(':')[1].split('-')[0]
				endBorder=infoCircu.split('|')[1].split(':')[1].split('-')[1].replace("+","")
				strand=infoCircu.split('|')[1].split(':')[1][-1]
				description=chromosome+'|'+str(startBorder)+'-'+str(endBorder)
				if (cGcC>=THRESHOLD_CGCC and g4H>=THRESHOLD_G4H and g4NN>=THRESHOLD_G4NN): # if window contain value > several threshold
					passed=True # passage over the differents threshold for this window
					if (oldPassed != passed): # if it's the first windows (beginning if the passage
						listeCGcC=[] # list which will contain the values of score cgCG for windows over the THRESHOLD_CGCC
						listeG4H=[]  # list which will contain the values of score G4H for windows over the THRESHOLD_G4H
						listeG4NN=[]  # list which will contain the values of score G4NN for windows over the THRESHOLD_G4NN
						descriptionOverThreshold=description # assignation of description for this variable
						sequenceG4=sequence
						oldPassed=passed # assignation, we pass over the thresholds
						if (strand == '+'): # if gene positif
							startG4=startWindow # startG4 same as startWindow of G4 screener
							endG4=endWindow # endG4 same as endWindow of G4 screener
						else:	# if gene negatif
							startG4=int(endBorder)-int(startWindow)+1 # calcul of startG4
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
					
						if (strand == '+'): # if gene positif
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
						if (strand == '+'):
							headerG4=idCircu+"|"+str(startG4)+"|"+str(endG4)+"|"+strand
						elif (strand == '-'):
							headerG4=idCircu+"|"+str(endG4)+"|"+str(startG4)+"|"+strand
					
						if (strand == '+'): # if gene positif
							if (int(startG4) > int(startBorder) and int(startG4) < int(endBorder) and int(endG4) > int(endBorder)):# overpass
								condition=True
						else:
							if (int(startG4) > int(startBorder) and int(startG4) < int(endBorder) and int(endG4) < int(startBorder)):# overpass
								condition=True
						if (condition==True):  # if overpass
							localisation='junction'
							if (strand == '+'): # if gene +
								endG4=int(startBorder) + (int(endG4) - int(endBorder)) # calcul new endG4
								headerG4=idCircu+"|"+str(startG4)+"|"+str(endG4)+"|"+strand
							else: # if gene -
								endG4=int(endBorder) - (int(startBorder) - int(endG4)) # calcul new endG4
								headerG4=idCircu+"|"+str(endG4)+"|"+str(startG4)+"|"+strand
							G4Detected[headerG4]=str(meanCGcC), str(meanG4Hunter),sequenceG4 , str(meanG4NN), localisation, rnaType
						else:
							if (G4Detected.has_key(headerG4) == False and strand != None):
								G4Detected[headerG4]=str(meanCGcC), str(meanG4Hunter),sequenceG4 , str(meanG4NN), localisation, rnaType
						condition=False
						localisation='ExonNC'
			
	return G4Detected


def ReturnLengthInCircu(Length, inputfile): ## for G4 in gene
	""" Add informations in dictionary of length of each circular regions

	    Parameters
	    ----------
	    Length
		dictionary of length of each circular regions
	    inputfile : string
		name of file which contain the liste of informations for each circular 
	    Returns
	    -------
	    Length
		dictionary of length of each circular regions
	"""
	inputfile= open(inputfile,"r")	# file opening for reading
	for line in inputfile:	# for each file csv create by G4 screener
		if (re.search('^[0-9]', line)): # if the line is not the header of the file 
			words=line.split('\t') # parsing by the separator, here '\t' (file csv)
			if (len(words) == 8):
				numRow=words[0].rstrip() #numero of row for this for this window
				infoCircu=words[1].rstrip()
				cGcC=float(words[2].rstrip()) #  score cGcC for this window
				g4H=float(words[3].rstrip()) #  score G4Hunter for this window 
				sequence=words[4].rstrip()	# sequence for this window
				startWindow=int(words[5].rstrip())	# start window from G4screener
				endWindow=int(words[6].rstrip()) ## end window from G4screener
				g4NN=float(words[7].rstrip()) # score G4NN for this window
				idCircu=infoCircu.split('|')[0] # id sequence circular
				chromosome=infoCircu.split('|')[1].split(':')[0] # chromosome
				startBorder=infoCircu.split('|')[1].split(':')[1].split('-')[0]
				endBorder=infoCircu.split('|')[1].split(':')[1].split('-')[1].replace("+","")
				strand=infoCircu.split('|')[1].split(':')[1][-1]
				description=chromosome+'|'+str(startBorder)+'-'+str(endBorder)
				if (Length.has_key(idCircu) == False ):
					Length[idCircu]=GetLengthFraction(startBorder,endBorder)
	return Length
		

######################################################################################################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	#parser.add_argument ('-i', '--inputfile', default = '/home/local/USHERBROOKE/bels2814/Documents/These/circulaire/Circu/human_hg19_circRNAs_putative_spliced_sequence_jonction_00001.csv')
	parser.add_argument ('-i', '--inputfile', default = '/home/local/USHERBROOKE/bels2814/Documents/These/circulaire/Circu/TEST.csv')
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


	G4DetectedInCircu={}
	LengthInCircu={}
	length=0
	
	parametersTool=[THRESHOLD_CGCC,THRESHOLD_G4H,THRESHOLD_G4NN, WINDOW, STEP]
		
	
	directory='/home/local/USHERBROOKE/bels2814/Documents/These/circulaire/Circu'
	for path, dirs, files in os.walk(directory): # for each element of the directory to passed
		for filename in files: # for each files
			inputfile=directory+'/'+filename
			if ('TEST' not in inputfile):
				print filename
				G4DetectedInCircu=ReturnG4InCircu(G4DetectedInCircu, inputfile, parametersTool) ### detection of PG4 in circular case
				LengthInCircu=ReturnLengthInCircu(LengthInCircu, inputfile)
	print '----------------------------------------'
	sommeLength=0
	for key, values in LengthInCircu.items():
		sommeLength=sommeLength+values
	print sommeLength , len(G4DetectedInCircu)
	print len(G4DetectedInCircu)/ float(sommeLength) * 1000


	print '---------------------------------------------'
	print '\nNumber PG4 detected in circRNA:	'+str(len(G4DetectedInCircu))
	print 'Number total circRNA tested:	'+str(len(LengthInCircu))
	print 'Length of all circRNA (nt):	'+str(sommeLength)
	print 'Densité in Kb :	'+round(len(G4DetectedInCircu)/ float(sommeLength) * 1000,2)
	print 'Percent circRNA with Kb:	'+str(round(len(G4DetectedInGene)/ float(len(piRNA)) * 100,2))


	
	
	
main()
	
