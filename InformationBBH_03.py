#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
import re
from pprint import pprint

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
	""" Create a dictionary with for each gene we will have the transcript id and
	the protein ID, but all this is only available if there is a protein ID 
	Parameters
	    ----------
	     filename : string, name of the file that contain for each gene the transcript id and the protein id
	    Returns
	    -------
	     DicoBBH : keys (string) -> rG4s from specie 1; value (string) -> rG4s from specie 2
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
	     DicoHomology : {"Gene":{idGeneSpecie1 : [idGeneSpecie2]}
						"Transcript" : {idTranscririptSpecie1 : [idTranscriptSpecie2]}}
	"""
	DicoHomology = {"Gene" : {}, "Transcript" : {}}

	Specie1GeneID, Specie1TranscriptID, Specie2GeneID, Specie2TranscriptID = GetIDTypeSpecie(dicoTypeIdSpecies, specie1, specie2)
			
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l :
				words=l.split('\t')
				GeneSpecie1 = words[0]
				TranscriptSpecie1 = words[1]
				GeneSpecie2 = words[2]
				TranscriptSpecie2 = words[3]
				
				if GeneSpecie1 not in DicoHomology["Gene"]: # if the gene from the specie 1 is not
					# in the dictionary, we add it and we create the list of homologs genes
					DicoHomology["Gene"][GeneSpecie1] = []
					DicoHomology["Gene"][GeneSpecie1].append(GeneSpecie2) # add of the "first" homologue gene
				elif GeneSpecie1 in DicoHomology["Gene"] and GeneSpecie2 not in DicoHomology["Gene"][GeneSpecie1] :
				 # the the gene from specie one is already in the dictionary but he 
				 # got a new homologue gene 
				 # (it is possible that the same pairs of homologue appear in the file, that's why this step is needed)
					DicoHomology["Gene"][GeneSpecie1].append(GeneSpecie2)
				
				# same but for transcript
				if re.search(Specie2TranscriptID, TranscriptSpecie2): # first we need to verify
					# if the homology is with a transcrit or a protein
					if TranscriptSpecie1 not in DicoHomology["Transcript"] :
						DicoHomology["Transcript"][TranscriptSpecie1] = [TranscriptSpecie2]
					elif GeneSpecie1 in DicoHomology["Transcript"] and TranscriptSpecie2 not in DicoHomology["Transcript"][TranscriptSpecie1] : 
						DicoHomology["Transcript"][TranscriptSpecie1].append(TranscriptSpecie2)
				else: # the id is the id of a protein
					if TranscriptSpecie2 in DicoGeneTranscriptProteinSpecie2["Protein"] :
						# if the prot id is in the dico GTP so get the id of the transcript
						TranscFromProtSpecie2 = DicoGeneTranscriptProteinSpecie2["Protein"][TranscriptSpecie2][0]
						if TranscriptSpecie1 not in DicoHomology["Transcript"] : 
							# if the transcript is not alreayd in the dictionary 
							# we add it with the first homologue encounter
							DicoHomology["Transcript"][TranscriptSpecie1] = [TranscFromProtSpecie2]
						elif GeneSpecie1 in DicoHomology["Transcript"] and TranscFromProtSpecie2 not in DicoHomology["Transcript"][TranscriptSpecie1] :
							# the transcript from the specie 1 is in the dico 
							# but we found an other homologue from the specie 2
							# so we add the last one to the liisst of homologue
							DicoHomology["Transcript"][TranscriptSpecie1].append(TranscFromProtSpecie2)
					#~ else : # the prot id is not in the dico GTP, that's not normal
					# bad anotation we don't keep them 

	return(DicoHomology)
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
	return(DicoInfoG4)
#----------------------------------------------------------------------#
def FindHomologyInGene(dicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2, DicoInfoG4):
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
		
	    Returns
	    -------
	    DicoInfoG4["Gene"] : {idBBH : {"Homology" : {idCoupleHomologue : {"Common Localisations" : [Localisations]
																			"Common Biotypes" : [Biotypes]}}}}
	"""
	for idG4Sp1 in dicoBBH :
		idG4Sp2 = dicoBBH[idG4Sp1]
		idBBH = str(idG4Sp1+"|"+idG4Sp2)
		
		# first we get the list of genes where the G4 of each specie is
		if idG4Sp1 in DicoInfoG4Sp1["Gene"] :
			genesSp1 =  DicoInfoG4Sp1["Gene"][idG4Sp1]
		else :
			print "Patate 2"
		if idG4Sp2 in DicoInfoG4Sp2["Gene"] :
			genesSp2 =  DicoInfoG4Sp2["Gene"][idG4Sp2]
		else :
			print "Patate 3"
		
		# Then we get homology
		for geneSp1 in genesSp1 :
			if geneSp1 in DicoHomology["Gene"] : # if this gene have some homologues
				HgeneSp1 = DicoHomology["Gene"][geneSp1]
			else : # there is no homologues for this gene
				HgeneSp1 = []
			
			if list(set(HgeneSp1).intersection(genesSp2)) : 
				# there an intersection between the homologues of the gene of specie 1
				# and with the genes of the specie 2 where the G4 is
				for geneSp2 in list(set(HgeneSp1).intersection(genesSp2)) :
					# for the first Specie
					if idG4Sp1 in DicoInfoG4Sp1["Gene"] :
						if geneSp1 in DicoInfoG4Sp1["Gene"][idG4Sp1] :
							localisationsSp1 = DicoInfoG4Sp1["Gene"][idG4Sp1][geneSp1]["Localisations"]
							biotypesSp1 = DicoInfoG4Sp1["Gene"][idG4Sp1][geneSp1]["Biotype"]
						else :
							print "Patate 4"
					else :
						print "Patate 5"
					
					# for the second Specie
					if idG4Sp2 in DicoInfoG4Sp2["Gene"] :
						if geneSp2 in DicoInfoG4Sp2["Gene"][idG4Sp2] :
							localisationsSp2 = DicoInfoG4Sp2["Gene"][idG4Sp2][geneSp2]["Localisations"]
							biotypesSp2 = DicoInfoG4Sp2["Gene"][idG4Sp2][geneSp2]["Biotype"]
						else :
							print "Patate 6"
					else :
						print "Patate 7"
					
					commonLocalisations = []
					commonBiotypes = []
					if list(set(localisationsSp1).intersection(localisationsSp2)) :
						commonLocalisations = list(set(localisationsSp1).intersection(localisationsSp2))
					if list(set(biotypesSp1).intersection(biotypesSp2)) :
						commonBiotypes = list(set(biotypesSp1).intersection(biotypesSp2))
					
					coupleHomologue = str(geneSp1+"|"+geneSp2)
					if idBBH not in DicoInfoBBH["Gene"] :
						# First time the couple of BBH is encounter
						# so we create all information relative to it 
						DicoInfoBBH["Gene"][idBBH] = {"Homology" : {}}
						DicoInfoBBH["Gene"][idBBH]["Homology"][coupleHomologue] = {}
						DicoInfoBBH["Gene"][idBBH]["Homology"][coupleHomologue]["Common Localisations"] = commonLocalisations
						DicoInfoBBH["Gene"][idBBH]["Homology"][coupleHomologue]["Common Biotypes"] = commonBiotypes
					elif idBBH in DicoInfoBBH["Gene"] and coupleHomologue not in DicoInfoBBH["Gene"][idBBH]["Homology"] :
						# the couple of BBH is already in the dictionary but 
						# not the couple there homologues's genes
						DicoInfoBBH["Gene"][idBBH]["Homology"][coupleHomologue] = {}
						DicoInfoBBH["Gene"][idBBH]["Homology"][coupleHomologue]["Common Localisations"] = commonLocalisations
						DicoInfoBBH["Gene"][idBBH]["Homology"][coupleHomologue]["Common Biotypes"] = commonBiotypes
	
	return(DicoInfoBBH["Gene"])
#----------------------------------------------------------------------#
def FindInfoHomology(dicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2):
	""" Create a dictonary for each G4 with his localisation
	Parameters
	    ----------
	     dicoBBH : dictionnary of BBH, see ImportBBH doc
	     DicoHomology : {"Gene":{idGeneSpecie1 : [idGeneSpecie2]}
						"Transcript" : {idTranscririptSpecie1 : [idTranscriptSpecie2]}}
	    Returns
	    -------
	"""
	DicoInfoBBH = {"Gene" : {}, "Transcript" : {}}
	
	
						
						
				
			
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
	#~ dicoGeneTranscriptProteinSpecie1 = GetIDGeneTranscriptProteins(filenameGeneTranscriptProteinSpecie2)
	#~ dicoGeneTranscriptProteinSpecie2 = GetIDGeneTranscriptProteins(filenameGeneTranscriptProteinSpecie2)
	#~ dicoHomology = ImportHomology(filenameHomology, dicoGeneTranscriptProteinSpecie2, dicoTypeIdSpecies, specie1, specie2)
	#~ dicoInfoG4Specie1 = importInfoG4(filenameInfoSpecie1, dicoGeneTranscriptProteinSpecie1)
	#~ dicoInfoG4Specie2 = importInfoG4(filenameInfoSpecie2, dicoGeneTranscriptProteinSpecie2)
	
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
	
	#~ dicoBBH, dicoHomology, dicoBiotypeSpecie1, dicoBiotypeSpecie2, dicoLocalisationSpecie1, dicoLocalisationSpecie2 = importData(path, specie1, specie2, idTypePerSpecie)
	importData(path, specie1, specie2, idTypePerSpecie)
	
	#~ infoBBH = FindCommonInformations(dicoBBHv, dicoLocalisationSpecie1, dicoLocalisationSpecie2, dicoBiotypeSpecie1, dicoBiotypeSpecie2)
	
	# write output specie 1 
	#~ output= open(path+specie2+'_BBHrG4.txt',"w") # file opening for reading
	#~ output.write("\n".join(BBHG4sSpecie2))
	#~ output.close()
#----------------------------------------------------------------------#	
main()
