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
			idProtein = words[2]
			if idGene not in DicoGeneTranscriptProtein :
				DicoGeneTranscriptProtein[idGene] = {idTranscrit : idProtein}
			else if idGene in DicoGeneTranscriptProtein and idTranscrit not in DicoGeneTranscriptProtein[idGene] :
				DicoGeneTranscriptProtein[idGene][idTranscrit] = idProtein
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
	DicoHomologyGene = {}
	DicoHomologyTranscript = {}

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
			
			if GeneSpecie1 not in DicoHomologyGene: # if the gene from the specie 1 is not
				# in the dictionary, we add it and we create the list of homologs genes
				DicoHomologyGene[GeneSpecie1] = []
				DicoHomologyGene[GeneSpecie1].append(GeneSpecie2) # add of the "first" homologue gene
			else if GeneSpecie1 in DicoHomologyGene and GeneSpecie2 not in DicoHomologyGene[GeneSpecie1]
			 # the the gene from specie one is already in the dictionary but he 
			 # got a new homologue gene 
			 # (it is possible that the same pairs of homologue appear in the file, that's why this step is needed)
				DicoHomologyGene[GeneSpecie1].append(GeneSpecie2)
			
			# same but for transcript
			if re.search(Specie2TranscriptID, TranscriptSpecie2): # first we need to verify
				# if the homology is with a transcrit or a protein
				if TranscriptSpecie1 not in DicoHomologyTranscript :
					DicoHomologyTranscript[TranscriptSpecie1] = []
					DicoHomologyTranscript[TranscriptSpecie1].append(TranscriptSpecie2)
				else if GeneSpecie1 in DicoHomologyTranscript and TranscriptSpecie2 not in DicoHomologyTranscript[TranscriptSpecie1]
					DicoHomologyTranscript[TranscriptSpecie1].append(TranscriptSpecie2)
			else:
				for i in DicoGeneTranscriptProteinSpecie2[GeneSpecie2]:
					if TranscriptSpecie2 == DicoGeneTranscriptProteinSpecie2[GeneSpecie2][i].value():
						idTranscriptSpecie2 = DicoGeneTranscriptProteinSpecie2[GeneSpecie2][i].value()
				if TranscriptSpecie1 not in DicoHomologyTranscript :
					DicoHomologyTranscript[TranscriptSpecie1] = []
					DicoHomologyTranscript[TranscriptSpecie1].append(TranscriptSpecie2)
				else if GeneSpecie1 in DicoHomologyTranscript and idTranscriptSpecie2 not in DicoHomologyTranscript[TranscriptSpecie1]
					DicoHomologyTranscript[TranscriptSpecie1].append(idTranscriptSpecie2)
	return(DicoHomologyGene, DicoHomologyTranscript)
#----------------------------------------------------------------------#
def importBiotype(filename):
	""" Create a dictonary for each gene, all his transcript and the biotype of the transcript
	Parameters
	    ----------
	     filename : string, name of the file containing all the gene and his transcript
	     with their biotype
	    Returns
	    -------
	     DicoBiotype : keys (string) -> geneID; value (dictionary) -> key = Transcript id, value = biotype
	"""
	DicoBiotype = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			words=line.split('|')
			gene = words[0]
			transcript = words[1]
			biotype = words[3]
			if gene not in DicoBiotype : # if the gene is no already in the dictionary we add it
				DicoBiotype[gene] = {transcript : biotype}
			else if gene in DicoBiotype and transcript not in DicoBiotype[gene] :
				#if the gene is already in the dictionary but not the transcript we add it 
				DicoBiotype[gene][transcript] = biotype
	return(DicoBiotype)
#----------------------------------------------------------------------#
def importLocalisation(filename):
	""" Create a dictonary for each G4 with his localisation
	Parameters
	    ----------
	     filename : string, name of the file containing for all the G4 in Genome their localisation
	    Returns
	    -------
	     DicoLocation : keys (string) -> geneID; value (list of string) -> all the localisation of a rG4s
	"""
	DicoLocalisation = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			words=line.split('\t')
			ids = words[0]
			ids = ids.split(':')
			G4 = str(ids[1]+":"+ids[2])
			G4 = G4.split('|')[0]
			localisations = words[1].split(';')
			if G4 not in DicoLocalisation : # if the gene is no already in the dictionary we add it
				DicoLocalisation[G4] = localisations
			else : # that's not suppose to happen 
				print "Patate 1"
	return(DicoLocalisation)
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
def FindCommonLocalisation(dicoBBH, DicoLocalisationSpecie1, DicoLocalisationSpecie2):
	""" import of all the data we need 
	Parameters
	    ----------
	    dicoBBH : dictionnary of BBH, see ImportBBH doc 
	    DicoLocalisationSpecie1 : dictionary of Localisation for specie 1, see importLocalisation
	    DicoLocalisationSpecie2 : dictionary of Localisation for specie 2, see importLocalisation
	    Returns
	    -------
	    dicoLocBBH : dictionary, keys = G4BBH (string) and value = list of localisation in common (string)
	"""
	dicoLocBBH = {}
	for G4Specie1 in dicoBBH : # browse of all the BBH
		G4Specie2 = dicoBBH[G4Specie1]
		if G4Specie1 in DicoLocalisationSpecie1 and G4Specie2 in DicoLocalisationSpecie2 :
			localisationG4Specie1 = DicoLocalisationSpecie1[G4Specie1]
			localisationG4Specie2 = DicoLocalisationSpecie2[G4Specie2]
			localisationCommon = list(set(localisationG4Specie1).intersection(localisationG4Specie2))
			BBH = str(G4Specie1+"|"+G4Specie2)
			dicoLocBBH[BBH] = localisationCommon
		else : # all the G4 chould be in the dico localisation
			#so we should never enter here
			print "Patate 2"
	return(dicoLocBBH)
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
	
	
						
	# write output specie 1 
	output= open(path+specie2+'_BBHrG4.txt',"w") # file opening for reading
	output.write("\n".join(BBHG4sSpecie2))
	output.close()
#----------------------------------------------------------------------#	
main()


	










