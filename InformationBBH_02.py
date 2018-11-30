#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
import re

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
			words=line.split('\t')
			DicoBBH[words[0]]= words[1]
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
	     DicoBBH : keys (string) -> rG4s from specie 1; value (string) -> rG4s from specie 2
	"""
	re1 = r+spece1+" gene"
	re2 = r+spece1+" transcript"
	re3 = r+spece2+" gene"
	re4 = r+spece2+" transcript"
	for i in dicoTypeIdSpecies :
		if re.search(re1, i) : # if key is specie 1 gene
			Specie1GeneID = dicoTypeIdSpecies[i]
		else if re.search(re2, i) : # if key is specie 1 gene
			Specie1TranscriptID = dicoTypeIdSpecies[i]
		else if re.search(re3, i) : # if key is specie 1 gene
			Specie2GeneID = dicoTypeIdSpecies[i]
		else if re.search(re4, i) : # if key is specie 1 gene
			Specie2TranscriptID = dicoTypeIdSpecies[i]
	return(Specie1GeneID, Specie1TranscriptID, Specie2GeneID, Specie2TranscriptID)
#----------------------------------------------------------------------#
def GetIDGeneTranscriptProteins(filename):
	""" Create a dictionary with for each gene we will have the transcript id and
	the protein ID, but all this is only available if there is a protein ID 
	Parameters
	    ----------
	     filename : string, name of the file that contain for each gene the transcript id and the protein id
	    Returns
	    -------
	     DicoBBH : keys (string) -> rG4s from specie 1; value (string) -> rG4s from specie 2
	"""
	DicoGeneTranscriptProtein = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			words=line.split('\t')
			idGene = words[0]
			idTranscrit = words[1]
			if words[2] : # there is a protein id 
				idProtein = words[2]
				if idTranscrit not in DicoGeneTranscriptProtein : # if the transcript has not been encounter yet
					# we add it to the dictionary with the first protein encounter
					DicoGeneTranscriptProtein[idTranscrit] = {idGene : [idProtein]}
				else if idTranscrit in DicoGeneTranscriptProtein and idProtein not in DicoGeneTranscriptProtein[idTranscrit][idGene] :
					# if the transcript is already in the dictionary but not the idprot
					DicoGeneTranscriptProtein[idTranscrit][idGene].append(idProtein)
			else : # there is no protein id
				idProtein = "No protein"
				if idTranscrit not in DicoGeneTranscriptProtein : # if the transcript has not been encounter yet
					# we add it to the dictionary with the first protein encounter
					DicoGeneTranscriptProtein[idTranscrit] = {idGene : [idProtein]}
				else if idTranscrit in DicoGeneTranscriptProtein and idProtein not in DicoGeneTranscriptProtein[idTranscrit][idGene] :
					# if the transcript is already in the dictionary but not the idprot
					DicoGeneTranscriptProtein[idTranscrit][idGene].append(idProtein)
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
	     DicoHomologyGene : keys (sring) -> Gene of specie 1
							value (list of string) -> Gene of specie 2 that 
								are homologue to the specie 1 for the gene in key
		 DicoHomologyTranscript : keys (sring) -> Transcript of specie 1
								 value (list of string) -> Transcript of specie 2 that 
								are homologue to the specie 1 for the gene in key
	"""
	DicoHomology = {"Gene" : {}, "Transcript" : {}}

	Specie1GeneID, Specie1TranscriptID, Specie2GeneID, Specie2TranscriptID = GetIDTypeSpecie(dicoTypeIdSpecies, specie1, specie2)
			
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			words=line.split('\t')
			GeneSpecie1 = words[0]
			TranscriptSpecie1 = words[1]
			GeneSpecie2 = words[2]
			TranscriptSpecie2 = words[3]
			
			if GeneSpecie1 not in DicoHomology["Gene"]: # if the gene from the specie 1 is not
				# in the dictionary, we add it and we create the list of homologs genes
				DicoHomology["Gene"][GeneSpecie1] = []
				DicoHomology["Gene"][GeneSpecie1].append(GeneSpecie2) # add of the "first" homologue gene
			else if GeneSpecie1 in DicoHomology["Gene"] and GeneSpecie2 not in DicoHomology["Gene"][GeneSpecie1]
			 # the the gene from specie one is already in the dictionary but he 
			 # got a new homologue gene 
			 # (it is possible that the same pairs of homologue appear in the file, that's why this step is needed)
				DicoHomology["Gene"][GeneSpecie1].append(GeneSpecie2)
			
			# same but for transcript
			if re.search(Specie2TranscriptID, TranscriptSpecie2): # first we need to verify
				# if the homology is with a transcrit or a protein
				if TranscriptSpecie1 not in DicoHomologyTranscript :
					DicoHomologyTranscript[TranscriptSpecie1] = []
					DicoHomologyTranscript[TranscriptSpecie1].append(TranscriptSpecie2)
				else if GeneSpecie1 in DicoHomology["Transcript"] and TranscriptSpecie2 not in DicoHomology["Transcript"][TranscriptSpecie1]
					DicoHomology["Transcript"][TranscriptSpecie1].append(TranscriptSpecie2)
			else: # the id is the id of a protein
				if TranscriptSpecie2 in DicoGeneTranscriptProteinSpecie2[GeneSpecie2] :
					# if the prot id is in the dico GTP so get the id of the transcript
					TranscFromProtSpecie2 = DicoGeneTranscriptProteinSpecie2[GeneSpecie2][TranscriptSpecie2]
					if TranscriptSpecie1 not in DicoHomology["Transcript"] : 
						# if the transcript is not alreayd in the dictionary 
						# we add it with the first homologue encounter
						DicoHomology["Transcript"][TranscriptSpecie1] = []
						DicoHomology["Transcript"][TranscriptSpecie1].append(TranscFromProtSpecie2)
					else if GeneSpecie1 in DicoHomology["Transcript"] and TranscFromProtSpecie2 not in DicoHomology["Transcript"][TranscriptSpecie1] :
						# the transcript from the specie 1 is in the dico 
						# but we found an other homologue from the specie 2
						# so we add the last one to the liisst of homologue
						DicoHomology["Transcript"][TranscriptSpecie1].append(TranscFromProtSpecie2)
				else : # the prot id is not in the dico GTP, that's not normal
					print "Patate 4"
	return(DicoHomologyGene, DicoHomologyTranscript)
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
	     DicoInfoG4 :	{"Transcript" : {idTranscript : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
						 "Gene" : {idGene : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
	"""
	DicoInfoG4 = {"Gene": {}, "Transcript" : {}}
	
	# from the file G4InTranscript we will extract the data of transcript : localisation and biotype
	with open(filenameG4InTranscript) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			words=line.split('\t')
			ids = words[0]
			ids = ids.split('|')
			idTranscript = ids[0]
			G4 = ids[1]
			localisation = words[5]
			biotype = words[6]
			if G4 not in DicoInfoG4["Transcript"]: # the g4 is not yet encounter
				DicoInfoG4["Transcript"][G4][idTranscript]["Localisations"] = [localisation]
				DicoInfoG4["Transcript"][G4][idTranscript]["Biotype"] = [biotype]
			else if G4 in DicoInfoG4["Transcript"] and idTranscript not in DicoInfoG4["Transcript"][G4] :
				# the G4 is encounter but not in this transcript
				DicoInfoG4["Transcript"][G4][idTranscript]["Localisations"] = [localisation]
				DicoInfoG4["Transcript"][G4][idTranscript]["Biotype"] = [biotype]
			else : # that's not suppose to happen 
				print "Patate 1"
	
	# with the new informations and the dicoGTP we can get the missing ones :
	# for genes level the localisation and biotypes
	
	for G4 in DicoInfoG4["Transcript"]: # browse all G4
		G4PerTranscript = DicoInfoG4["Transcript"][G4] # get the sub dictionary of transcript where this G4 is
		for Transcrit in G4PerTranscript : # browse all the transcript of this G4
			gene = DicoGTP[Transcrit].keys() # get only the name of the transcript's gene
			if G4 not in DicoInfoG4["Gene"] : # if the G4 is not already in the dictionary
				# we add it with the first gene encounter and his relative informations
				DicoInfoG4["Gene"][G4][gene] = G4PerTranscript[Transcrit]
			else if G4 in DicoInfoG4["Gene"] and gene not in DicoInfoG4["Gene"][G4] :
				# if the G4 is already in the dictionary but not the gene 
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
			else if G4 in DicoInfoG4["Gene"] and gene in DicoInfoG4["Gene"][G4] :
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
	return(DicoInfoG4)
#----------------------------------------------------------------------#
def FindCommonInformationsGene(dicoBBH, DicoInfoG4Specie1, DicoInfoG4Specie2):
	""" find common informations between two BBH at the level of genome
	Parameters
	    ----------
	    dicoBBH : dictionnary of BBH, see ImportBBH doc 
	    DicoInfoG4Specie1 : dictionary of G4's information  for specie 1, see importInfoG4
	    DicoInfoG4Specie2 : dictionary of G4's information  for specie 2, see importInfoG4
	    Returns
	    -------
	    dicoInformationCommon : 
	"""
	dicoInformationCommon = {"Location" : {}, "Biotype" : {}}
	for G4Specie1 in dicoBBH : # browse of all the BBH
		G4Specie2 = dicoBBH[G4Specie1]
		
		if G4Specie1 in DicoInfoG4Specie1["Gene"] and G4Specie2 in DicoInfoG4Specie2["Gene"] :
			infoG4Specie1 = DicoInfoG4Specie1["Gene"][G4]
			infoG4Specie2 = DicoInfoG4Specie2["Gene"][G4]
			localisationCommon = list(set(localisationG4Specie1).intersection(localisationG4Specie2))
			BBH = str(G4Specie1+"|"+G4Specie2)
			dicoInformationCommon["Location"][BBH] = localisationCommon

	return(dicoInformationCommon)
#----------------------------------------------------------------------#
def FindInfoHomology(dicoBBH, DicoHomologyGene, DicoHomologyTranscript):
	""" Create a dictonary for each G4 with his localisation
	Parameters
	    ----------
	     dicoBBH : dictionnary of BBH, see ImportBBH doc
	     DicoHomologyGene : keys (sring) -> Gene of specie 1
							value (list of string) -> Gene of specie 2 that 
								are homologue to the specie 1 for the gene in key
		 DicoHomologyTranscript : keys (sring) -> Transcript of specie 1
								 value (list of string) -> Transcript of specie 2 that 
								are homologue to the specie 1 for the gene in key
	    Returns
	    -------
	"""
	DicoHomologyCommon = {"Gene" : {}, "Transcript" : {}}
	for G4Specie1 in dicoBBH :
		G4Specie2 = dicoBBH[G4Specie1]
		
		# First at the gene level
		if 
			
	return()
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
	     	dicoBiotypeSpecie1
	     	dicoBiotypeSpecie2
	     	dicoLocalisationSpecie1
	     	dicoLocalisationSpecie2

	"""
	filenameBBH = path+'BBH.txt'
	filenameGeneTranscriptProteinSpecie1 = path+specie1+"_All_GTP.txt"
	filenameHomoly = path+specie1+"_homoloie_"+specie2+"_All.txt"
	filenameBiotypeSpecie1 = path+specie1+"_All_TranscriptType.txt"
	filenameLocalisationSpecie1 = path+specie1+"_All_G4InGenome.txt"
	filenameBiotypeSpecie2 = path+specie2+"_All_TranscriptType.txt"
	filenameLocalisationSpecie2 = path+specie2+"_All_G4InGenome.txt"
	
	dicoBBH = ImportBBH(filenameBBH)
	dicoGeneTranscriptProteinSpecie2 = GetIDGeneTranscriptProteins(filenameGeneTranscriptProteinSpecie2)
	dicoHomology = ImportHomology(filenameHomoly, dicoGeneTranscriptProteinSpecie2)
	dicoBiotypeSpecie1 = importBiotype(filenameBiotypeSpecie1)
	dicoBiotypeSpecie2 = importBiotype(filenameBiotypeSpecie2)
	dicoLocalisationSpecie1 = importLocalisation(filenameLocalisationSpecie1)
	dicoLocalisationSpecie2 = importLocalisation(filenameLocalisationSpecie2)
	
	return(dicoBBH, dicoHomology, dicoBiotypeSpecie1, dicoBiotypeSpecie2, dicoLocalisationSpecie1, dicoLocalisationSpecie2)
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
	
	dicoBBH, dicoHomology, dicoBiotypeSpecie1, dicoBiotypeSpecie2, dicoLocalisationSpecie1, dicoLocalisationSpecie2 = importData(path, specie1, specie2, idTypePerSpecie)
	
	infoBBH = FindCommonInformations(dicoBBHv, dicoLocalisationSpecie1, dicoLocalisationSpecie2, dicoBiotypeSpecie1, dicoBiotypeSpecie2)
	
	# write output specie 1 
	output= open(path+specie2+'_BBHrG4.txt',"w") # file opening for reading
	output.write("\n".join(BBHG4sSpecie2))
	output.close()
#----------------------------------------------------------------------#	
main()
