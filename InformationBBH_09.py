#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
import re
import ImportBBH as imp
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
def AddInformationsOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleOrtho, commonLocalisations, commonBiotypes, typeInfo):
	""" Adds the informations to the final dictionary containing All 
		infromation about the BBH if there is some orthology 
		Parameters
	    ----------
	    DicoInfoBBHLevel : {idBBH : {typeInfo : {couple : {"Common Localisations" : [commonLocalisations]
														   "Common Biotypes" : commonBiotypes}}}
		idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		coupleOrtho : string, id of the couple like geneSp1|geneSp2 or transcriptSp1|transcriptSp2
		commonLocalisations : list of string, contain the common localisation 
							  between the gene/transcripts that are in the list coupleOrtho.
							  This list could be empty
		commonBiotypes : list of string, contain the common biotypes 
						 between the gene/transcripts that are in the list coupleOrtho.
						 This list could be empty
		typeInfo : string, here it is "Orthology"
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {"Orthology" : {coupleOrtho : {"Common Localisations" : [commonLocalisations]
																	   "Common Biotypes" : commonBiotypes}}}}
	"""
	idBBH = str(idG4Sp1+"|"+idG4Sp2)
	if idBBH not in DicoInfoBBHLevel :
		# First time the couple of BBH is encounter
		# so we create all information relative to it 
		DicoInfoBBHLevel[idBBH] = {typeInfo : {}}
		DicoInfoBBHLevel[idBBH][coupleOrtho] = {}
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and coupleOrtho not in DicoInfoBBHLevel[idBBH][typeInfo] :
		# the couple of BBH is already in the dictionary but 
		# not the couple there homologues's levelS
		DicoInfoBBHLevel[idBBH][coupleOrtho] = {}
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and coupleOrtho in DicoInfoBBHLevel[idBBH] :
		# the couple of BBH is already in the dictionary and also the couple of homologue
		# yet we still have a chance to get some new informations so we need to check it
		if list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)) :
			# if there is new localisations we add them
			DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Localisations"].append(list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)))
		if list(set(DicoInfoBBHLevel[idBBH]["Common Biotypes"]) ^ set(commonBiotypes)) :
			# if there is new biotypes we add them
			DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Biotypes"].append(commonBiotypes)
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def AddInfoNoOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, newList, commonLocalisations, commonBiotypes, typeInfo):
	""" Adds the informations to the final dictionary containing All infromation 
		about the BBH if there is no orthology links
		Parameters
	    ----------
	    DicoInfoBBHLevel : {idBBH : {typeInfo : {couple : {"Common Localisations" : [commonLocalisations]
														   "Common Biotypes" : commonBiotypes}}}
		idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		newList : list of string, id of the gene that are not orthologues
						but they still have common informations in common
		commonLocalisations : list of string, contain the common localisation 
							  between the gene/transcripts that are in the list coupleOrtho.
							  This list could be empty
		commonBiotypes : list of string, contain the common biotypes 
						 between the gene/transcripts that are in the list coupleOrtho.
						 This list could be empty
		typeInfo : string, here it is "No orthology"
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {"Orthology" : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
																	   "Common Biotypes" : commonBiotypes}}}}
	"""
	idBBH = str(idG4Sp1+"|"+idG4Sp2)
	newList = tuple(newList)
	if idBBH not in DicoInfoBBHLevel :
		# First time the couple of BBH is encounter
		# so we create all information relative to it
		DicoInfoBBHLevel[idBBH] = {typeInfo : {}}
		DicoInfoBBHLevel[idBBH][typeInfo][newList] = {}
		DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and newList not in DicoInfoBBHLevel[idBBH][typeInfo] :
		cpt = 0
		# the couple of BBH is already in the dictionary but 
		# not the id that are not orthologues, this can be du to many reasons
		# maybe there just one id in plus the list or maybe the list doesn't exist
		for existingList in DicoInfoBBHLevel[idBBH][typeInfo]:
			if set(existingList).issubset(set(newList)) :
				cpt +=1
				# the existingList is a subset of newList
				newComLoca = list(set(commonLocalisations) ^ (DicoInfoBBHLevel[idBBH][typeInfo][existingList]["Common Localisations"]))
				newComBiot = list(set(commonBiotypes) ^ (DicoInfoBBHLevel[idBBH][typeInfo][existingList]["Common Biotypes"]))
				DicoInfoBBHLevel[idBBH][typeInfo][newList] = dictionary.pop(DicoInfoBBHLevel[idBBH][typeInfo][existingList])
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"].append(newComLoca)
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"].append(newComBiot)
			else :
				# the list of id is new
				DicoInfoBBHLevel[idBBH][typeInfo][newList] = {}
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"] = commonLocalisations
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"] = commonBiotypes
		if cpt > 1 :
			print "Patate 7"		
	elif idBBH in DicoInfoBBHLevel and newList in DicoInfoBBHLevel[idBBH][typeInfo] :
		# the couple of BBH is already in the dictionary and also the couple of homologue
		# yet we still have a chance to get some new informations so we need to check it
		if list(set(DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"]) ^ set(commonLocalisations)) :
			# if there is new localisations we add them
			DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"].append(list(set(DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"]) ^ set(commonLocalisations)))
		if list(set(DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"]) ^ set(commonBiotypes)) :
			# if there is new biotypes we add them
			DicoInfoBBHLevel[idBBH][idNoOrthology][typeInfo][newList]["Common Biotypes"].append(list(set(DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"]) ^ set(commonBiotypes)))
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetInfoNoOrthology(levelSp1, nonHomologueslevelSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel):
	""" Here we test if there is common information between gene that hasn't orthology links
		Parameters
	    ----------
	    nonHomologueslevelSp2 : list of string, list of the gene that are not orthologs to the gene of 
								the specie 1 where there is the G4 BBH we are looking at for
		levelSp1 : string, transcript/gene from the specie 1
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleNoOrth : {"Common Localisations" : [commonLocalisations]
													 "Common Biotypes" : commonBiotypes}}}
	"""
	idNoOrthology = [levelSp1]
	commonLocalisationsNoOrth = []
	commonBiotypesNoOrth[]
	for levelSp2 in nonHomologueslevelSp2 : 
		# retrivial of the localisations and biotype of the two specie
		# from the gene/transcript that are homologues
		localisationsSp1WithoutH, biotypesSp1WithoutH = GetLocalisationsAndBiotype(idG4Sp1, levelSp1, DicoInfoG4Sp1Level)
		localisationsSp2WithoutH, biotypesSp2WithoutH = GetLocalisationsAndBiotype(idG4Sp2, levelSp2, DicoInfoG4Sp2Level)
		commonLocalisationsWithoutH, commonBiotypesWithoutH = GetCommonInformation(localisationsSp1WithoutH, localisationsSp1WithoutH, biotypesSp1WithoutH, biotypesSp1WithoutH)
		if commonLocalisationsWithoutH or commonBiotypesWithoutH :
			commonLevelSp2.append(levelSp2)
			commonLocalisationsNoOrth.append(commonLocalisationsWithoutH)
			commonBiotypesNoOrth.append(commonBiotypesWithoutH)
	commonLocalisationsNoOrth = list(set(commonLocalisationsNoOrth))
	commonBiotypesNoOrth =  list(set(commonBiotypesNoOrth))
	DicoInfoBBHLevel = AddInfoNoOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, newList, commonLocalisations, commonBiotypes, "No orthology")
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetInfoOfParalogues(listParalogues, idG4s,DicoHomologyLevel, DicoInfoG4Level):
	""" The gene are orthologues, so we will search if they have also paralogues
		Parameters
	    ----------
	    listParalogues : list of string, list that contain all the paralogues of a gene/transcript
	    idG4 : string, id of the G4s that we are trying to analise
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
	    Returns
	    -------
	    paraloguesG4 : {"Localisations" : {Loc1 : [paralogues]}
						"Biotypes" : {Bt1 : [paralogues]}}
	"""

	paraloguesG4 : {"Localisations" : {}, "Biotypes" : {}}
	for paralogue in listParalogues :
		if paralogue in DicoInfoG4Level[idG4s] :
			for loc in DicoInfoG4Level[idG4s]["Localisations"]:
				if loc not in paraloguesG4["Localisations"] :
					paraloguesG4["Localisations"][loc] = []
					paraloguesG4["Localisations"][loc].append(paralogue)
				elif loc in paraloguesG4["Localisations"] and paralogue not in paraloguesG4["Localisations"][loc] : 
					paraloguesG4["Localisations"][loc].append(paralogue)
			for Bt in DicoInfoG4Level[idG4s]["Biotypes"]
				if Bt not in paraloguesG4["Biotypes"] :
					paraloguesG4["Biotypes"][Bt] = []
					paraloguesG4["Biotypes"][Bt].append(paralogue)
				elif Bt in paraloguesG4["Biotypes"] and paralogue not in paraloguesG4["Biotypes"][Bt] : 
					paraloguesG4["Biotypes"][Bt].append(paralogue)
		else :
			# the paralogue of the gene/transcript don't have the G4s
			print str(paralogue) + " est un paralogue mais n'a pas le pG4s de sont paralogue."
	return(paraloguesG4)
#----------------------------------------------------------------------#
def GetCommonInfoOfParalogues(paraloguesG4Sp1, paraloguesG4Sp2):
	""" The gene are orthologues, so we will search if they have also paralogues
		Parameters
	    ----------
	    paraloguesG4Sp1/Sp2 : {"Localisations" : {Loc1 : [paralogues]}
							   "Biotypes" : {Bt1 : [paralogues]}}
	    Returns
	    -------
	    commonParalogues : list of pralogues that have common information
	    commonloc : list of common localisations
	    commonBt : list of common biotype
	"""
	commonloc = []
	commonBt = []
	commonParalogues = []
	for loc in paraloguesG4Sp1["Localisations"] :
		if loc in paraloguesG4Sp2["Localisations"] :
			# the localisation is present in both species's paralogues
			# so we add those information to the common lists
			commonParalogues.append(paraloguesG4Sp2["Localisations"][loc])
			commonParalogues.append(paraloguesG4Sp1["Localisations"][loc])
			commonloc.append(loc)
	for Bt in paraloguesG4Sp1["Biotypes"] :
		if Bt in paraloguesG4Sp2["Biotypes"] :
			# the biotype is present in both species's paralogues
			# so we add those information to the common lists
			commonParalogues.append(paraloguesG4Sp2["Biotypes"][Bt])
			commonParalogues.append(paraloguesG4Sp1["Biotypes"][Bt])
			commonBt.append(Bt)
	return(commonParalogues, commonloc, commonBt)
#----------------------------------------------------------------------#
def GetInfoParalogyWithOrtho(idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevelHom, DicoHomologyLevel):
	""" The gene are orthologues, so we will search if they have also paralogues
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
		
	    Returns
	    -------
	    DicoHomologyLevelHom : {commonParalogues : {"Common Localisations" : [commonLocalisations]
													"Common Biotypes" : [commonBiotypes]}}}
	"""
	if levelSp2 in DicoHomologyLevel]["Paralogues Sp2"] :
		paraloguesSp2 = DicoHomologyLevel]["Paralogues Sp2"][levelSp2]
		paraloguesSp2.append(levelSp2)
	if levelSp1 in DicoHomologyLevel]["Paralogues Sp1"] :
		paraloguesSp1 = DicoHomologyLevel]["Paralogues Sp1"][levelSp1]
		paraloguesSp1.append(levelSp1)
		
	paraloguesG4Sp1 = GetInfoOfParalogues(paraloguesSp1, idG4Sp1,DicoHomologyLevel, DicoInfoG4Sp1Level)
	paraloguesG4Sp2 = GetInfoOfParalogues(paraloguesSp2, idG4Sp2,DicoHomologyLevel, DicoInfoG4Sp2Level)
	
	commonParalogues, commonloc, commonBt = GetCommonInfoOfParalogues(paraloguesG4Sp1, paraloguesG4Sp2)
	
	if tuple(commonParalogues) not in DicoInfoBBHLevelHom :
		DicoInfoBBHLevelHom[tuple(commonParalogues)] = {"Common localisations" : [], "Common biotypes" : []}
		DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common localisations"].append(commonloc)
		DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common biotypes"].append(commonBt)
	else if tuple(commonParalogues) in DicoInfoBBHLevelHom :
		if set(commonloc) ^ DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common localisations"] :
			DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common localisations"].append(set(commonloc) ^ DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common localisations"])
		if set(commonBt) ^ DicoHomologyLevelHom[tuple(commonParalogues)]["Common biotypes"] :
			DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common biotypes"].append(set(commonBt) ^ DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common biotypes"])
	return(DicoInfoBBHLevelHom)
#----------------------------------------------------------------------#
def GetInfoHomology(HlevelSp1, levelSp1, levelsSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, DicoHomologyLevel):
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
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
	"""
	if list(set(HlevelSp1).intersection(levelsSp2)) :
		# if there an intersection between the homologues of the gene/transcript of specie 1
		# and with the levelS of the specie 2 where the G4 is
		
		for levelSp2 in list(set(HlevelSp1).intersection(levelsSp2)) :
			# retrivial of the localisations and biotype of the two specie
			# from the gene/transcript that are homologues
			localisationsSp1WithH, biotypesSp1WithH = GetLocalisationsAndBiotype(idG4Sp1, levelSp1, DicoInfoG4Sp1Level)
			localisationsSp2WithH, biotypesSp2WithH = GetLocalisationsAndBiotype(idG4Sp2, levelSp2, DicoInfoG4Sp2Level)
			commonLocalisationsWithH, commonBiotypesWithH = GetCommonInformation(localisationsSp1WithH, localisationsSp2WithH, biotypesSp1WithH, biotypesSp2WithH)
			coupleOrthologue = str(levelSp1+"|"+levelSp2)
			idBBH = str(idG4Sp1+"|"+idG4Sp2)
			DicoInfoBBHLevel.update(AddInformationsOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleOrthologue, commonLocalisationsWithH, commonBiotypesWithH, "Orthology"))
			if idBBH not in DicoInfoBBHLevel[idBBH] :
				DicoInfoBBHLevel[idBBH] = {"Orthology" : {}, "Homology" : {}}
			DicoInfoBBHLevel[idBBH]["Homology"].update(GetInfoParalogyWithOrtho(idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel[idBBH]["Homology"], DicoHomologyLevel))
			
		nonHomologueslevelSp2 = list(set(HlevelSp1) ^ set(levelsSp2))
		if nonHomologueslevelSp2 : 
		# retrivial of the list of G4 that are not homologues (unique in the two lists)
			DicoInfoBBHLevel = GetInfoNoOrthology(levelSp1, nonHomologueslevelSp2, idG4Sp2, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel)
			
	else :
		# there is no element in comon => no orthologues between the gene/transcript
		# where the G4 BBH are
		DicoInfoBBHLevel = GetInfoNoOrthology(levelSp1, levelsSp2, idG4Sp2, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel)
return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetDicoInfoBBH(idG4Sp1, idG4Sp2, DicoHomologyLevel, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, level)
	""" Create a dictonary for each G4 with his localisation and biotype,
		the informations that are retrieve depend on the level transcript/gene
		Parameters
	    ----------
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
							 "Paralogy Sp1" : {idGeneSp1 : [idGeneSp1]}
							 "Paralogy Sp2" : {idGeneSp2 : [idGeneSp2]}}
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
	"""
	# first we retrieve the dico of informations from a level (transcript/gene)
	# where the G4 of each specie is
	if idG4Sp1 in DicoInfoG4Sp1Level : # retrievial of informations for the specie 1
		levelsSp1 =  DicoInfoG4Sp1Level[idG4Sp1]
	else :
		print "Patate 2 : G4 from specie 1 doesn't exist in the dictionary of G4"
		
	if idG4Sp2 in DicoInfoG4Sp2Level : # retrievial of informations for the specie 2
		levelsSp2 =  DicoInfoG4Sp2Level[idG4Sp2]
	else :
		print "Patate 3 : G4 from specie 2 doesn't exist in the dictionary of G4"
			
	# Then we browse all genes where the G4 from the specie 1 is to get
	# all the informations we want (Paralogy  Sp1/Sp2, Orthology or not)
	for levelSp1 in levelsSp1 : # browse all genes where the G4 from the specie 1 is
		DicoOrthologyLevel = DicoHomologyLevel["Orthologues"][level]
		DicoParalogyLevel = DicoHomologyLevel["Paralogues"][level]
		
		# Orthology
		if levelSp1 in DicoOrthologyLevel : # if this gene have some orthologues
			OrtLevelSp1 = DicoOrthologyLevel[levelSp1] # retrevial of the 
			# gene/transcript from the specie 2 that are orthologues to levelSp1
		else : # there is no orthologues for this gene
			OrtLevelSp1 = []
		DicoInfoBBHLevel = GetInfoHomology(OrtLevelSp1, levelSp1, levelsSp2, idG4Sp2, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, DicoHomologyLevel)
		
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetDicoInfo(dicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2, DicoInfoBBHLevel, level):
	""" Create a dictonary for each G4 with his localisation and biotype
		the informations that are retrieve depend on the level transcript/gene
	Parameters
	    ----------
	     dicoBBH : dictionnary of BBH, see ImportBBH doc
	     DicoHomology : {"Gene": {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								  "Paralogy Sp1" : {idGeneSp1 : [idGeneSp1]}
								  "Paralogy Sp2" : {idGeneSp2 : [idGeneSp2]}}
						"Transcript" : {"Orthology" : {idTranscririptSp1 : [idTranscriptSp2]}}
										"Paralogy Sp1" : {{idTranscririptSp1 : [idTranscriptSp1]}}
										"Paralogy Sp2" : {{idTranscririptSp2 : [idTranscriptSp2]}} 
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
	
	dicoBBH = imp.ImportBBH(filenameBBH)
	print "Import of BBH -> Done"
	dicoGeneTranscriptProteinSpecie1 = imp.ImportIDGeneTranscriptProteins(filenameGeneTranscriptProteinSpecie2)
	print "Import of GTP from specie 1 -> Done"
	dicoGeneTranscriptProteinSpecie2 = imp.ImportIDGeneTranscriptProteins(filenameGeneTranscriptProteinSpecie2)
	print "Import of GTP from specie 2 -> Done"
	dicoHomology = imp.ImportHomology(filenameHomology, dicoGeneTranscriptProteinSpecie2, dicoTypeIdSpecies, specie1, specie2)
	print "Import of Homology -> Done"
	dicoInfoG4Specie1 = imp.importInfoG4(filenameInfoSpecie1, dicoGeneTranscriptProteinSpecie1)
	print "Import of informations from specie 1 -> Done"
	dicoInfoG4Specie2 = imp.importInfoG4(filenameInfoSpecie2, dicoGeneTranscriptProteinSpecie2)
	print "Import of informations from specie 1 -> Done"
	
	#~ pprint(dicoGeneTranscriptProteinSpecie2)
	#~ pprint(dicoHomology)
	#~ pprint(dicoInfoG4Specie2)
	return(dicoBBH, dicoHomology, dicoGeneTranscriptProteinSpecie2, dicoInfoG4Specie1, dicoInfoG4Specie2)
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
