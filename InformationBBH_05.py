#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
import re
from pprint import pprint

#----------------------------------------------------------------------#
def PrintCounter(number, length):
	""" Print where the computation of the principal dictionary (DicoInfoBBH) is
	Parameters
	    ----------
	     number : the number of BBH passed
	     length : number total of BBH
	"""
	bool10 = False
    bool20 = False
    bool30 = False
    bool40 = False
    bool50 = False
    bool60 = False
    bool70 = False
    bool80 = False
    bool90 = False
    bool100 = False
    stat = round((float(number)/float(length))*100,1)
	if stat == 10 and bool10 == False:
		bool10 = True
		print "-10%"
	elif stat == 20 and bool20 == False:
		bool20 = True
		print "--20%"
	elif stat == 30 and bool30 == False:
		bool30 = True
		print "---30%"
	elif stat == 40 and bool40 == False:
		bool40 = True
		print "----40%"
	elif stat == 50 and bool50 == False:
		bool50 = True
		print "-----50%"
	elif stat == 60 and bool60 == False:
		bool60 = True
		print "------60%"
	elif stat == 70 and bool70 == False:
		bool70 = True
		print "-------70%"
	elif stat == 80 and bool80 == False:
		bool80 = True
		print "--------80%"
	elif stat == 90 and bool90 == False:
		bool90 = True
		print "---------90%"
	elif stat == 100 and bool100 == False:
		bool100 = True
		print "----------100%"
#----------------------------------------------------------------------#
def ImportBBH(filename):
	""" Create a dictonary that contain all the regions G4s and there BBH
	Parameters
	    ----------
	     filename : string, name of the file containing all the gene
	     of both specie as BBH
	    Returns
	    -------
	     DicoBBH : keys (string) -> rG4s from specie 1; value (string) -> rG4s from specie 2
	"""
	DicoBBH = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			words=l.split('\t')
			idG4Sp1 = words[0] # the first word correspond to the id of the G4 from the first specie
			idG4Sp2 = words[1] # the second word correspond to the G4 from the second specie
			DicoBBH[idG4Sp1]= idG4Sp2  # we add the two id of the G4 in a dictionary
	return(DicoBBH)
#----------------------------------------------------------------------#
def GetIDTypeSpecie(dicoTypeIdSpecies, specie1, specie2):
	""" find what is the type of the ensemble ID for transcript and gene for our two specie
	Parameters
	    ----------
	     dicoTypeIdSpecies : keys (string) -> description of the id type like HS gene
							 values (string) -> id type like ENSG
	     specie1 : string, name of the first specie like HS for Homo sapiens
	     specie2 : string, name of the first specie like MM for Mus musculus
	    Returns
	    -------
	     Specie1GeneID : string, type of the id of the gene for the specie 1 (ex : ENSMUSG00000038564)
	     Specie1TranscriptID : string, type of the id of the transcript for the specie 1 (ex : ENSMUST00000063414)
	     Specie2GeneID : string, type of the id of the gene for the specie 2 (ex : ENSMUSG00000038564)
	     Specie2TranscriptID : string, type of the id of the transcript for the specie 2 (ex : ENSMUST00000063414)
	"""
	re1 = specie1+" gene"
	re2 = specie1+" transcript"
	re3 = specie2+" gene"
	re4 = specie2+" transcript"
	for i in dicoTypeIdSpecies :
		if re.search(re1, i) : # if key is specie 1 gene
			Specie1GeneID = dicoTypeIdSpecies[i]
		elif re.search(re2, i) : # if key is specie 1 gene
			Specie1TranscriptID = dicoTypeIdSpecies[i]
		elif re.search(re3, i) : # if key is specie 1 gene
			Specie2GeneID = dicoTypeIdSpecies[i]
		elif re.search(re4, i) : # if key is specie 1 gene
			Specie2TranscriptID = dicoTypeIdSpecies[i]
	return(Specie1GeneID, Specie1TranscriptID, Specie2GeneID, Specie2TranscriptID)
#----------------------------------------------------------------------#
def GetIDGeneTranscriptProteins(filename):
	""" Create a dictionary with for each gene we will have the transcript ID and
	the protein ID, if there is no protein ID in ENSEMBL there will be a label "No protein"
	Parameters
	    ----------
	     filename : string, name of the file that contain for each gene the transcript id and the protein id
	    Returns
	    -------
	     DicoGeneTranscriptProtein : {"Transcript" : {idTranscript : {idGene : idProtein}}
									  "Protein" : {idProtein : idTranscript}}
	"""
	DicoGeneTranscriptProtein = {"Transcript" : {}, "Protein" : {}}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l : 
				words=l.split('\t')
				#~ print l
				idGene = words[0]
				idTranscrit = words[1]
				if words[2] : # there is a protein id 
					idProtein = words[2]
					if idTranscrit not in DicoGeneTranscriptProtein["Transcript"] : # if the transcript has not been encounter yet
						# we add it to the dictionary with the first protein encounter
						DicoGeneTranscriptProtein["Transcript"][idTranscrit] = {idGene : [idProtein]}
					elif idTranscrit in DicoGeneTranscriptProtein["Transcript"] and idProtein not in DicoGeneTranscriptProtein["Transcript"][idTranscrit][idGene] :
						# if the transcript is already in the dictionary but not the idprot
						DicoGeneTranscriptProtein["Transcript"][idTranscrit][idGene].append(idProtein)
					if idProtein not in DicoGeneTranscriptProtein["Protein"] : # if the transcript has not been encounter yet
						# we add it to the dictionary with the first protein encounter
						DicoGeneTranscriptProtein["Protein"][idProtein] = [idTranscrit]
					elif idProtein in DicoGeneTranscriptProtein["Protein"] and idTranscrit not in DicoGeneTranscriptProtein["Protein"][idProtein] :
						# if the transcript is already in the dictionary but not the idprot
						DicoGeneTranscriptProtein["Protein"][idProtein].append(idTranscrit)
				else : # there is no protein id
					idProtein = "No protein"
					if idTranscrit not in DicoGeneTranscriptProtein["Transcript"] : # if the transcript has not been encounter yet
						# we add it to the dictionary with the first protein encounter
						DicoGeneTranscriptProtein["Transcript"][idTranscrit] = {idGene : [idProtein]}
					elif idTranscrit in DicoGeneTranscriptProtein["Transcript"] and idProtein not in DicoGeneTranscriptProtein["Transcript"][idTranscrit][idGene] :
						# if the transcript is already in the dictionary but not the idprot
						DicoGeneTranscriptProtein["Transcript"][idTranscrit][idGene].append(idProtein)
					if idProtein not in DicoGeneTranscriptProtein["Protein"] : # if the transcript has not been encounter yet
						# we add it to the dictionary with the first protein encounter
						DicoGeneTranscriptProtein["Protein"][idProtein] = [idTranscrit]
					elif idProtein in DicoGeneTranscriptProtein["Protein"] and idTranscrit not in DicoGeneTranscriptProtein["Protein"][idProtein] :
						# if the transcript is already in the dictionary but not the idprot
						DicoGeneTranscriptProtein["Protein"][idProtein].append(idTranscrit)
	return(DicoGeneTranscriptProtein)
#----------------------------------------------------------------------#
def ImportHomology(filename, DicoGeneTranscriptProteinSpecie2, dicoTypeIdSpecies, specie1, specie2):
	""" Create two dictonary that contain Homologie between the two specie,
	one with the homology between genes of the two specie and the other one
	between the transcripts of the two specie
	Parameters
	    ----------
	     filename : string, name of the file containing all the gene and transcript
	     that are homologue between two specie
	    Returns
	    -------
	     DicoHomology : {"Gene": {"Paralogues" : {idGeneSpecie1 : [idGeneSpecie1]}
								  "Orthologues" : {idGeneSpecie1 : [idGeneSpecie2]}}
						"Transcript" : {"Paralogues" : {idTranscririptSpecie1 : [idTranscriptSpecie1]}
										"Orthologues" : {idTranscririptSpecie1 : [idTranscriptSpecie2]}}
	"""
	DicoHomology = {"Gene" : {"Paralogues" : {}, "Orthologues" : {}}, "Transcript" : {"Paralogues" : {}, "Orthologues" : {}}}

	Specie1GeneID, Specie1TranscriptID, Specie2GeneID, Specie2TranscriptID = GetIDTypeSpecie(dicoTypeIdSpecies, specie1, specie2)
			
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l :
				words=l.split('\t')
				GeneSpecie1 = words[0] # gene of the specie 1
				TranscriptSpecie1 = words[1] # transcript/protein of the specie 1
				GeneSpecie1Paralogue = words[2] # gene id of a paralogue of GeneSpecie1
				TranscriptSpecie1Paralogue = words[3] # transcript id of a paralogue of TranscriptSpecie1
				GeneSpecie2 = words[4] # orthologue gene from the specie 2
				TranscriptSpecie2 = words[5] # orthologue transcript from the specie 2
				
				if GeneSpecie1 not in DicoHomology["Gene"]["Orthologues"]: # if the gene from the specie 1 is not
					# in the dictionary, we add it and we create the list of homologs genes
					DicoHomology["Gene"]["Orthologues"][GeneSpecie1] = []
					DicoHomology["Gene"]["Orthologues"][GeneSpecie1].append(GeneSpecie2) # add of the "first" homologue gene
				elif GeneSpecie1 in DicoHomology["Gene"]["Orthologues"] and GeneSpecie2 not in DicoHomology["Gene"]["Orthologues"][GeneSpecie1] :
				 # the the gene from specie one is already in the dictionary but he 
				 # got a new homologue gene 
				 # (it is possible that the same pairs of homologue appear in the file, that's why this step is needed)
					DicoHomology["Gene"]["Orthologues"][GeneSpecie1].append(GeneSpecie2)
				
				# same but for transcript
				if re.search(Specie2TranscriptID, TranscriptSpecie2): # first we need to verify
					# if the homology is with a transcrit or a protein
					if TranscriptSpecie1 not in DicoHomology["Transcript"]["Orthologues"] :
						DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] = [TranscriptSpecie2]
					elif GeneSpecie1 in DicoHomology["Transcript"]["Orthologues"] and TranscriptSpecie2 not in DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] : 
						DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1].append(TranscriptSpecie2)
				else: # the id is the id of a protein
					if TranscriptSpecie2 in DicoGeneTranscriptProteinSpecie2["Protein"] :
						# if the prot id is in the dico GTP so get the id of the transcript
						TranscFromProtSpecie2 = DicoGeneTranscriptProteinSpecie2["Protein"][TranscriptSpecie2][0]
						if TranscriptSpecie1 not in DicoHomology["Transcript"]["Orthologues"] : 
							# if the transcript is not alreayd in the dictionary 
							# we add it with the first homologue encounter
							DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] = [TranscFromProtSpecie2]
						elif GeneSpecie1 in DicoHomology["Transcript"]["Orthologues"] and TranscFromProtSpecie2 not in DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] :
							# the transcript from the specie 1 is in the dico 
							# but we found an other homologue from the specie 2
							# so we add the last one to the liisst of homologue
							DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1].append(TranscFromProtSpecie2)
					#~ else : # the prot id is not in the dico GTP, that's not normal
					# bad anotation we don't keep them 

	return(DicoHomology)
#----------------------------------------------------------------------#
def importInfoG4Transcript(filenameG4InTranscript):
	""" Create a dictonary for each G4 with his localisations and biotype at the 
		transcript level
	Parameters
	    ----------
	     filenameG4InTranscript : string, name of the file containing all 
	     the G4 In Transcript with their localisations and biotype (transcript level)
	    Returns
	    -------
	     DicoInfoG4["Transcript"] :	{idG4 : {idTranscript : {"Localisation" : [localisation]
															 "Biotype" : [Biotype]}}}
	"""
	# from the file G4InTranscript we will extract the data of transcript : localisation and biotype
	with open(filenameG4InTranscript) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l :
				words=l.split('\t')
				ids = words[0]
				ids = ids.split('|')
				idTranscript = ids[0]
				G4 = ids[1]
				localisation = words[5]
				biotype = words[6]
				if G4 not in DicoInfoG4["Transcript"]: # the g4 is not yet encounter
					DicoInfoG4["Transcript"][G4] = {}
					DicoInfoG4["Transcript"][G4][idTranscript] = {"Localisations" : []}
					DicoInfoG4["Transcript"][G4][idTranscript]["Localisations"].append(localisation)
					DicoInfoG4["Transcript"][G4][idTranscript]["Biotype"] = [biotype]
				elif G4 in DicoInfoG4["Transcript"] and idTranscript not in DicoInfoG4["Transcript"][G4] :
					# the G4 is encounter but not in this transcript
					DicoInfoG4["Transcript"][G4][idTranscript] = {"Localisations" : []}
					DicoInfoG4["Transcript"][G4][idTranscript]["Localisations"] = [localisation]
					DicoInfoG4["Transcript"][G4][idTranscript]["Biotype"] = [biotype]
				else : # that's not suppose to happen 
					print "Patate 1"
	return(DicoInfoG4["Transcript"])
#----------------------------------------------------------------------#
def importInfoG4Gene(DicoGTP):
	""" Create a dictonary for each G4 with his localisations and biotype at the 
		gene level
	Parameters
	    ----------
	     DicoGTP : {"Transcript" : {idTranscript : {idGene : idProtein}}
									  "Protein" : {idProtein : idTranscript}}
	    Returns
	    -------
	     DicoInfoG4 : {idG4 : {idGene : {"Localisation" : [localisation]
										 "Biotype" : [Biotype]}}}
	"""
	# with the new informations and the dicoGTP we can get the missing ones :
	# for genes level the localisation and biotypes
	
	for G4 in DicoInfoG4["Transcript"]: # browse all G4
		G4PerTranscript = DicoInfoG4["Transcript"][G4] # get the sub dictionary of transcript where this G4 is
		for Transcrit in G4PerTranscript : # browse all the transcript of this G4
			gene = DicoGTP["Transcript"][Transcrit].keys() # get only the name of the transcript's gene
			if len(gene) > 1 :
				print "PATATE"
			gene = gene[0]
			if G4 not in DicoInfoG4["Gene"] : # if the G4 is not already in the dictionary
				# we add it with the first gene encounter and his relative informations
				DicoInfoG4["Gene"][G4] = {}
				DicoInfoG4["Gene"][G4][gene] = {}
				DicoInfoG4["Gene"][G4][gene] = G4PerTranscript[Transcrit]
			elif G4 in DicoInfoG4["Gene"] and gene not in DicoInfoG4["Gene"][G4] :
				# if the G4 is already in the dictionary but not the gene
				DicoInfoG4["Gene"][G4][gene] = {}
				DicoInfoG4["Gene"][G4][gene] = G4PerTranscript[Transcrit]
			elif G4 in DicoInfoG4["Gene"] and gene in DicoInfoG4["Gene"][G4] :
				# if the G4 is already in the dictionary but and the gene also
				# but for some other transcript we may have other informations 
				# we add it after verifying if we have new informations
				localisationGene = DicoInfoG4["Gene"][G4][gene]["Localisations"]
				biotypeGene = DicoInfoG4["Gene"][G4][gene]["Biotype"]
				localisationTranscript = G4PerTranscript[Transcrit]["Localisations"]
				biotypeTranscript = G4PerTranscript[Transcrit]["Biotype"]
				# 2 if because there can those two condition can be verified
				if localisationTranscript not in localisationGene :
					DicoInfoG4["Gene"][G4][gene]["Localisations"].append(localisationTranscript)
				if biotypeTranscript not in biotypeGene : 
					DicoInfoG4["Gene"][G4][gene]["Biotype"].append(biotypeTranscript)
	return(DicoInfoG4["Gene"])
#----------------------------------------------------------------------#
def importInfoG4(filenameG4InTranscript, DicoGTP):
	""" Create a dictonary for each G4 with his localisation
	Parameters
	    ----------
	     filenameG4InTranscript : string, name of the file containing all 
	     the G4 In Transcript with their localisations and biotype (transcript level)
	     DicoGTP : keys (string) -> rG4s from specie 1; value (string) -> rG4s from specie 2
	    Returns
	    -------
	     DicoInfoG4 :	{"Transcript" : {idG4 : {idTranscript : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
						 "Gene" : {idG4 : {idGene : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
	"""
	DicoInfoG4 = {"Gene": {}, "Transcript" : {}}
	
	DicoInfoG4["Transcript"] = importInfoG4Transcript(filenameG4InTranscript)
	DicoInfoG4["Gene"] = importInfoG4Gene(DicoGTP)
	
	return(DicoInfoG4)
#----------------------------------------------------------------------#
def AddingInformations(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleHomologue, commonLocalisations, commonBiotypes):
	""" Adds the informations to the final dictionary containing All infromation about the BBH
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
	idBBH = str(idG4Sp1+"|"+idG4Sp2) # creation of the id of the couple of BBH
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
		# the couple of BBH is already in the dictionary and also the couple of homologue
		# yet we still have a chance to get some new informations so we 
		# need to check it
		if list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)) :
			# if there is new localisations we add them
			DicoInfoBBHLevel[idBBH][coupleHomologue]["Common Localisations"].append(list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)))
		if list(set(DicoInfoBBHLevel[idBBH]["Common Biotypes"]) ^ set(commonBiotypes)) :
			# if there is new biotypes we add them
			DicoInfoBBHLevel[idBBH][coupleHomologue]["Common Biotypes"].append(commonBiotypes)
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetCommonInformation(localisationsSp1, localisationsSp2, biotypesSp1, biotypesSp2):
	""" Get the commons informations between two list for the biotypes and the localisations,
		if there is no commons informations, the list will stay empty
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
	""" Here we test many conditions (with or without homology), then we 
		retrieve the informations relative to them
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
	""" Create a dictonary for each G4 with his localisation and biotype
		the informations that are retrieve depend on the level transcript/gene
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
def GetDicoInfo(dicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2, DicoInfoBBHLevel, level):
	""" Create a dictonary for each G4 with his localisation and biotype
		the informations that are retrieve depend on the level transcript/gene
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
	cpt = 0
	lentot = len(dicoBBH)
	
	for idG4Sp1 in dicoBBH : # browse all BBH couple
		idG4Sp2 = dicoBBH[idG4Sp1]
		
		DicoInfoBBHLevel = GetDicoInfoBBH(idG4Sp1, idG4Sp2, DicoHomologyLevel, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, level)
		
		cpt += 1
		
		PrintCounter(cpt, len(dicoBBH)) # will print wheree  we are in the computation of this function
	
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def WriteInfoHomology(DicoInfoBBH, path, level):
	output = open(path+'Homology_'+level+'_BBHrG4.txt',"w") # file opening for reading
	output.write("BBH\tcouple Homologue\tLocalisations\tBiotypes")
	for BBH in DicoInfoBBH[level] :
		for coupleHomologue in DicoInfoBBH[level][BBH] :
			if coupleHomologue != "No Homology" :
				output.write(str(BBH)+"\t"+str(coupleHomologue)+"\t"+"|".join(DicoInfoBBH[level][BBH][coupleHomologue]["Common Localisations"])+"\t"+"|".join(DicoInfoBBH[level][BBH][coupleHomologue]["Common Biotypes"]))
	output.close()
	
#----------------------------------------------------------------------#
def WriteInfoWithoutHomology(DicoInfoBBH, path, level):
	output = open(path+'NoHomology_'+level+'_BBHrG4.txt',"w") # file opening for reading
	output.write("BBH\tcouple Homologue\tLocalisations\tBiotypes")
	for BBH in DicoInfoBBH[level] :
		for coupleHomologue in DicoInfoBBH[level][BBH] :
			output.write(str(BBH)+"\t"+"|".join(DicoInfoBBH[level][BBH]["No Homology"]["Common Localisations"])+"\t"+"|".join(DicoInfoBBH[level][BBH]["No Homology"]["Common Biotypes"]))
	output.close()
#----------------------------------------------------------------------#
def importData(path, specie1, specie2, dicoTypeIdSpecies):
	""" import of all the data we need 
	Parameters
	    ----------
	     filename : strin, name of the path for all file
	    Returns
	    -------
	     	dicoBBH
	     	dicoHomology
	     	dicoGeneTranscriptProteinSpecie2
	     	dicoInfoG4Specie1
	     	dicoInfoG4Specie2

	"""
	filenameBBH = path+'BBH.txt'
	filenameGeneTranscriptProteinSpecie1 = path+specie1+"_All_GTP.txt"
	filenameGeneTranscriptProteinSpecie2 = path+specie2+"_All_GTP.txt"
	filenameHomology = path+specie1+"_homologie_"+specie2+"_All.txt"
	filenameInfoSpecie1 = path+specie1+"_All_G4InTranscript.txt"
	filenameInfoSpecie2 = path+specie2+"_All_G4InTranscript.txt"
	
	#~ dicoBBH = ImportBBH(filenameBBH)
	# ~ print "Import of BBH -> Done"
	#~ dicoGeneTranscriptProteinSpecie1 = GetIDGeneTranscriptProteins(filenameGeneTranscriptProteinSpecie2)
	# ~ print "Import of GTP from specie 1 -> Done"
	#~ dicoGeneTranscriptProteinSpecie2 = GetIDGeneTranscriptProteins(filenameGeneTranscriptProteinSpecie2)
	# ~ print "Import of GTP from specie 2 -> Done"
	#~ dicoHomology = ImportHomology(filenameHomology, dicoGeneTranscriptProteinSpecie2, dicoTypeIdSpecies, specie1, specie2)
	# ~ print "Import of Homology -> Done"
	#~ dicoInfoG4Specie1 = importInfoG4(filenameInfoSpecie1, dicoGeneTranscriptProteinSpecie1)
	# ~ print "Import of informations from specie 1 -> Done"
	#~ dicoInfoG4Specie2 = importInfoG4(filenameInfoSpecie2, dicoGeneTranscriptProteinSpecie2)
	# ~ print "Import of informations from specie 1 -> Done"
	
	#~ pprint(dicoGeneTranscriptProteinSpecie2)
	#~ pprint(dicoHomology)
	#~ pprint(dicoInfoG4Specie2)
	#~ return(dicoBBH, dicoHomology, dicoGeneTranscriptProteinSpecie2, dicoInfoG4Specie1, dicoInfoG4Specie2)
#----------------------------------------------------------------------#
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'InfoBBH')
	parser.add_argument ('-p', '--path', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/Homologie/BBH/')
	parser.add_argument ('-sp1', '--specie1', default = 'HS')
	parser.add_argument ('-sp2', '--specie2', default = 'MM')
	return parser
#----------------------------------------------------------------------#
def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path	# directory which contain all the file with BBH and rG4s
	specie1=arg.specie1	# first specie to analyse
	specie2=arg.specie2	# first specie to analyse
	
	idTypePerSpecie = {"HS gene" : "ENSG0", "HS transcript" : "ENST0", "MM gene" : "ENSMUSG0", "MM transcript" : "ENSMUST0"}
	
	#~ dicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2 = importData(path, specie1, specie2, idTypePerSpecie)
	importData(path, specie1, specie2, idTypePerSpecie)
	
	#~ infoBBH = FindCommonInformations(dicoBBH, dicoLocalisationSpecie1, dicoLocalisationSpecie2, dicoBiotypeSpecie1, dicoBiotypeSpecie2)
	
	DicoInfoBBH = {"Gene" : {}, "Transcript" : {}}
	DicoInfoBBH["Gene"] = GetDicoInfo(dicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2, DicoInfoBBH, "Gene")
	DicoInfoBBH["Transcript"] = GetDicoInfo(dicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2, DicoInfoBBH, "Transcript")
	
	WriteInfoHomology(DicoInfoBBH, path, "Gene")
	print "File with homology for the gene level done"
	WriteInfoHomology(DicoInfoBBH, path, "Transcript")
	print "File with homology for the transcript level done"
	WriteInfoWithoutHomology(DicoInfoBBH, path, "Gene")
	print "File without homology for the gene level done"
	WriteInfoWithoutHomology(DicoInfoBBH, path, "Transcript")
	print "File without homology for the transcript level done"
#----------------------------------------------------------------------#	

main()
