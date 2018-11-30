#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
import re
from pprint import pprint

#----------------------------------------------------------------------#
def AddingInformations(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleHomologue, commonLocalisations, commonBiotypes):
	""" add the informations to the final dictionary containing All infromation about the BBH
		Parameters
	    ----------
	    DicoInfoBBHLevel :
		idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		coupleHomologue : 
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
							Yet, the coupleHomologue must be equal to "No homologues"
	"""
	idBBH = str(idG4Sp1+"|"+idG4Sp2)
	if idBBH not in DicoInfoBBHLevel :
		# First time the couple of BBH is encounter
		# so we create all information relative to it 
		DicoInfoBBHLevel[idBBH] = {typeInfo : {}}
		DicoInfoBBHLevel[idBBH][coupleHomologue] = {}
		DicoInfoBBHLevel[idBBH][coupleHomologue]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][coupleHomologue]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and coupleHomologue not in DicoInfoBBHLevel[idBBH][typeInfo] :
		# the couple of BBH is already in the dictionary but 
		# not the couple there homologues's levelS
		DicoInfoBBHLevel[idBBH][coupleHomologue] = {}
		DicoInfoBBHLevel[idBBH][coupleHomologue]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][coupleHomologue]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and coupleHomologue in DicoInfoBBHLevel[idBBH] :
		if list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)) :
			DicoInfoBBHLevel[idBBH][coupleHomologue]["Common Localisations"].append(list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)))
		if list(set(DicoInfoBBHLevel[idBBH]["Common Biotypes"]) ^ set(commonBiotypes)) :
			DicoInfoBBHLevel[idBBH][coupleHomologue]["Common Biotypes"].append(commonBiotypes)
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetCommonInformation(localisationsSp1, localisationsSp2, biotypesSp1, biotypesSp2):
	""" Get the commons informations between two list for the biotypes and the localisations
		Parameters
	    ----------
	    localisationsSp1 :	list of localisations from the specie 1
		localisationsSp2 :	list of localisations from the specie 2
		biotypesSp1 : list of biotypes from the specie 1
		biotypesSp2 : list of biotypes from the specie 2
		
	    Returns
	    -------
	    commonLocalisations : list of commons localisations of a G4 in a gene/transcript between the 2 specie
	    commonBiotypes : list of commons biotypes of a G4 in a gene/transcript between the 2 specie
	"""
	commonLocalisations = []
	commonBiotypes = []
	if list(set(localisationsSp1).intersection(localisationsSp2)) :
		commonLocalisations = list(set(localisationsSp1).intersection(localisationsSp2))
	if list(set(biotypesSp1).intersection(biotypesSp2)) :
		commonBiotypes = list(set(biotypesSp1).intersection(biotypesSp2))
	return(commonLocalisations, commonBiotypes)
#----------------------------------------------------------------------#
def GetLocalisationsAndBiotype(element1, element2, Dico):
	""" Get the localisations and the biotypes in a dictionary of dictionary 
		thanks to 2 elements
		Parameters
	    ----------
	    element1 :	string, element that bust be present in the first level of the dictionary
		element2 :	string, element that bust be present in the second level of the dictionary
		Dico :	{idTranscript/Gene : {idG4 : {"Localisation" : [localisation]
											  "Biotype" : [Biotype]}}}
	    Returns
	    -------
	    localisations : list of localisations of a G4 in a gene/transcript
	    biotypes : list of biotypes of a G4 in a gene/transcript
	"""
	if element1 in Dico :
		if element2 in Dico[element1] :
			localisations = Dico[element1][element2]["Localisations"]
			biotypes = Dico[element1][element2]["Biotype"]
		else :
			print "Patate 4"
	else :
		print "Patate 5"
	return(localisations, biotypes)
#----------------------------------------------------------------------#
def GetInfoHomology(HlevelSp1, levelSp1, levelsSp2, idG4Sp2, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel):
	""" Create a dictonary for each G4 with his localisation
		Parameters
	    ----------
	    HlevelSp1 : list of string, list of the gene that are homologue to the gene of 
					the specie 1 where there is the G4 BBH we are looking at for
		levelSp1 : string, transcript/gene from the specie 1
		levelsSp2 : list of string, list of transcript/gene from the specie 2
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		level : string, to choose if we want the gene level or the transcript level
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
							Yet, the coupleHomologue must be equal to "No homologues"
	"""
	if list(set(HlevelSp1).intersection(levelsSp2)) :
		# if there an intersection between the homologues of the gene of specie 1
		# and with the levelS of the specie 2 where the G4 is
		
		for levelSp2 in list(set(HlevelSp1).intersection(levelsSp2)) :
			# retrivial of the localisations and biotype of the two specie
			# from the gene/transcript that are homologues
			localisationsSp1WithH, biotypesSp1WithH = GetLocalisationsAndBiotype(idG4Sp1, levelSp1, DicoInfoG4Sp1Level)
			localisationsSp2WithH, biotypesSp2WithH = GetLocalisationsAndBiotype(idG4Sp2, levelSp2, DicoInfoG4Sp2Level)
			commonLocalisationsWithH, commonBiotypesWithH = GetCommonInformation(localisationsSp1WithH, localisationsSp2WithH, biotypesSp1WithH, biotypesSp2WithH)
			coupleHomologue = str(levelSp1+"|"+levelSp2)
			DicoInfoBBHLevel = AddingInformations(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleHomologue, commonLocalisationsWithH, commonBiotypesWithH)
			
		nonHomologueslevelSp2 = list(set(HlevelSp1) ^ set(levelsSp2)) 
		# retrivial of the list of G4 that are not homologues (unique in the two lists)
		for levelSp2 in nonHomologueslevelSp2 : 
			# retrivial of the localisations and biotype of the two specie
			# from the gene/transcript that are homologues
			localisationsSp1WithoutH, biotypesSp1WithoutH = GetLocalisationsAndBiotype(idG4Sp1, levelSp1, DicoInfoG4Sp1Level)
			localisationsSp2WithoutH, biotypesSp2WithoutH = GetLocalisationsAndBiotype(idG4Sp2, levelSp2, DicoInfoG4Sp2Level)
			commonLocalisationsWithoutH, commonBiotypesWithoutH = GetCommonInformation(localisationsSp1WithoutH, localisationsSp1WithoutH, biotypesSp1WithoutH, biotypesSp1WithoutH)
			coupleHomologue = "No Homology"
			DicoInfoBBHLevel = AddingInformations(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleHomologue, commonLocalisationsWithoutH, commonBiotypesWithoutH)
			
	else :
		# there is no element in comon => no homologues between the gene/transcript
		# where the G4 BBH are
		for levelSp2 in levelsSp2 :
			localisationsSp1WithoutH, biotypesSp1WithoutH = GetLocalisationsAndBiotype(idG4Sp1, levelSp1, DicoInfoG4Sp1Level)
			localisationsSp2WithoutH, biotypesSp2WithoutH = GetLocalisationsAndBiotype(idG4Sp2, levelSp2, DicoInfoG4Sp2Level)
			commonLocalisationsWithoutH, commonBiotypesWithoutH = GetCommonInformation(localisationsSp1WithoutH, localisationsSp1WithoutH, biotypesSp1WithoutH, biotypesSp1WithoutH)
			coupleHomologue = "No Homology"
			DicoInfoBBHLevel = AddingInformations(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleHomologue, commonLocalisationsWithoutH, commonBiotypesWithoutH)
return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetDicoInfoBBH(idG4Sp1, idG4Sp2, DicoHomologyLevel, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, level)
	""" Create a dictonary for each G4 with his localisation
		Parameters
	    ----------
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		level : string, to choose if we want the gene level or the transcript level
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
							Yet, the coupleHomologue must be equal to "No homologues"
	"""
	# first we get the dico of informations from a level (transcript/gene)
	# where the G4 of each specie is
	if idG4Sp1 in DicoInfoG4Sp1Level :
		levelsSp1 =  DicoInfoG4Sp1Level[idG4Sp1]
	else :
		print "Patate 2"
	if idG4Sp2 in DicoInfoG4Sp2Level :
		levelsSp2 =  DicoInfoG4Sp2Level[idG4Sp2]
	else :
		print "Patate 3"
			
	# Then we get common information with homology or without homology
	for levelSp1 in levelsSp1 : # browse all genes where the G4 from the specie 1 is
		if levelSp1 in DicoHomologyLevel : # if this gene have some homologues
			HlevelSp1 = DicoHomologyLevel[levelSp1]
		else : # there is no homologues for this gene
			HlevelSp1 = []
		
		DicoInfoBBHLevel = GetInfoHomology(HlevelSp1, levelSp1, levelsSp2, idG4Sp2, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel)
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetDicoInfo(dicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2, DicoInfoBBH, level):
	""" Create a dictonary for each G4 with his localisation
	Parameters
	    ----------
	     dicoBBH : dictionnary of BBH, see ImportBBH doc
	     DicoHomology : {"Gene":{idGeneSpecie1 : [idGeneSpecie2]}
						"Transcript" : {idTranscririptSpecie1 : [idTranscriptSpecie2]}}
		 DicoInfoG4Sp1 :{"Transcript" : {idTranscript : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
						 "Gene" : {idGene : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
		 DicoInfoG4Sp2 :{"Transcript" : {idTranscript : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
						 "Gene" : {idGene : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		level : string, to choose if we want the gene level or the transcript level
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
							Yet, the coupleHomologue must be equal to "No homologues"
	"""
	DicoInfoG4Sp1Level = DicoInfoG4Sp1[level]
	DicoInfoG4Sp2Level = DicoInfoG4Sp2[level]
	DicoHomologyLevel = DicoHomology[level]
	DicoInfoBBHLevel = DicoInfoBBH[level]
	
	for idG4Sp1 in dicoBBH : # browse all BBH couple
		idG4Sp2 = dicoBBH[idG4Sp1]
		
		DicoInfoBBHLevel = GetDicoInfoBBH(idG4Sp1, idG4Sp2, DicoHomologyLevel, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, level)
	
	return(DicoInfoBBH)
