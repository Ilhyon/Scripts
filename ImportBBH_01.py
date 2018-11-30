#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

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
def ImportIDGeneTranscriptProteins(filename):
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
def GetElement(words, numberElem):
	#~ print len(words), numberElem
	if len(words)-1 >= numberElem :
		element = words[numberElem]
	else :
		element = ""
	return(element)
#----------------------------------------------------------------------#
def InitialisationHomology(line, currentSp):
	""" Get the different element that we have in the file, it depend on
		which specie is the file
	Parameters
	    ----------
	     line : string, one line of the inputfile
	     currentSp : string, correspond to the current sp of the file
	    Returns
	    -------
	     GeneSpecie1 : string, gene of the specie 1
	     TranscriptSpecie1 : string, transcript/protein of the specie 1
	     GeneSpecie2Paralogue : string, gene id of a paralogue of GeneSpecie2
	     TranscriptSpecie2Paralogue : string, transcript id of a paralogue of TranscriptSpecie2
	     GeneSpecie1Paralogue : string, gene id of a paralogue of GeneSpecie1
	     TranscriptSpecie1Paralogue : string, transcript id of a paralogue of TranscriptSpecie1
	     GeneSpecie2 : string, orthologue gene from the specie 2
	     TranscriptSpecie2 : string, orthologue transcript from the specie 2
	"""
	words=line.split('\t')
	if currentSp == "Sp1" :
		GeneSpecie1 = GetElement(words, 0) # gene of the specie 1
		TranscriptSpecie1 = GetElement(words, 1) # transcript/protein of the specie 1
		GeneSpecie1Paralogue = GetElement(words, 2) # gene id of a paralogue of GeneSpecie1
		TranscriptSpecie1Paralogue = GetElement(words, 3) # transcript id of a paralogue of TranscriptSpecie1
		GeneSpecie2 = GetElement(words, 4) # orthologue gene from the specie 2
		TranscriptSpecie2 = GetElement(words, 5) # orthologue transcript from the specie 2
		return(GeneSpecie1, TranscriptSpecie1, GeneSpecie1Paralogue, TranscriptSpecie1Paralogue, GeneSpecie2, TranscriptSpecie2)
	else:
		GeneSpecie1 = GetElement(words, 4) # gene of the specie 1
		TranscriptSpecie1 = GetElement(words, 5) # transcript/protein of the specie 1
		GeneSpecie2Paralogue = GetElement(words, 2) # gene id of a paralogue of GeneSpecie2
		TranscriptSpecie2Paralogue = GetElement(words, 3) # transcript id of a paralogue of TranscriptSpecie2
		GeneSpecie2 = GetElement(words, 0) # orthologue gene from the specie 2
		TranscriptSpecie2 = GetElement(words, 1) # orthologue transcript from the specie 2
		return(GeneSpecie1, TranscriptSpecie1, GeneSpecie2Paralogue, TranscriptSpecie2Paralogue, GeneSpecie2, TranscriptSpecie2)
#----------------------------------------------------------------------#
def MajOrthology(GeneSpecie1, TranscriptSpecie1, GeneSpecie2, TranscriptSpecie2, DicoHomology):
	""" When the DicoHomology is already created and fill with the information 
		of the specie 1, we check if there is some update with the file of the specie 2
	Parameters
	    ----------
	     GeneSpecie1 : string, gene of the specie 1
	     TranscriptSpecie1 : string, transcript/protein of the specie 1
	     GeneSpecie2 : string, orthologue gene from the specie 2
	     TranscriptSpecie2 : string, orthologue transcript from the specie 2
	     DicoHomology : 
	    Returns
	    -------
	     DicoHomology : {"Gene": {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								  "Paralogy Sp1" : {idGeneSp1 : [idGeneSp1]}
								  "Paralogy Sp2" : {idGeneSp2 : [idGeneSp2]}}
						"Transcript" : {"Orthology" : {idTranscririptSp1 : [idTranscriptSp2]}}
										"Paralogy Sp1" : {{idTranscririptSp1 : [idTranscriptSp1]}}
										"Paralogy Sp2" : {{idTranscririptSp2 : [idTranscriptSp2]}} 
	"""
	
	# we start with the gene level
	if GeneSpecie1 in DicoHomology["Gene"]["Orthologues"] :
		if GeneSpecie2 not in DicoHomology["Gene"]["Orthologues"][GeneSpecie1] :
			#the gene from the specie 1 is known but there is a new orthologues form him
			DicoHomology["Gene"]["Orthologues"][GeneSpecie1].append(GeneSpecie2)
	else :
		#the gene from the specie 1 is not known so we add it
		DicoHomology["Gene"]["Orthologues"][GeneSpecie1] = []
		DicoHomology["Gene"]["Orthologues"][GeneSpecie1].append(GeneSpecie2)
	
	# then the transcript level
	if TranscriptSpecie1 in DicoHomology["Transcript"]["Orthologues"] :
		if TranscriptSpecie2 not in DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] :
			#the transcript from the specie 1 is known but there is a new orthologues form him
			DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1].append(TranscriptSpecie2)
	else :
		#the transcript from the specie 1 is not known so we add it
		DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] = []
		DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1].append(TranscriptSpecie2)
	return(DicoHomology)
#----------------------------------------------------------------------#
def AddGeneHomology(elem1, elem2, DicoHomologyLevel, homologyType):
	""" Add to the dictionary of homology (DicoHomology) the relation of orthology/paralogy,
		if those relations are not in the dictionary or even if they are 
		already in the dictionary, new information about the Gene/transcript are added
	Parameters
	    ----------
	     elem1 : string, correspond to a gene of the specie 1
	     elem2 : string, correspond to a gene of the specie 1 (for paralogy)
				 or to a gene of the specie 2 (for orthology)
	     DicoHomologyLevel : the subdictionary of DicoHomology with only 
							the level of gene or transcript
	     homologyType : string, can be "Paralogues" or "Orthologues"
	    Returns
	    -------
	     DicoHomologyLevel[homologyType][elem1] : {homologyType: {elem1 : [elem2]}}
	"""
	if elem1 not in DicoHomologyLevel: # if the gene from the specie 1 is not
		# in the dictionary, we add it and we create the list of homologs genes
		DicoHomologyLevel[homologyType][elem1] = []
		DicoHomologyLevel[homologyType][elem1].append(elem2) # add of the "first" homologue gene
	elif elem1 in DicoHomologyLevel[homologyType] and elem2 not in DicoHomologyLevel[homologyType][elem1] :
	 # the the gene from specie one is already in the dictionary but he 
	 # got a new homologue gene 
	 # (it is possible that the same pairs of gene homologue appear in the file, that's why this step is needed)
		DicoHomologyLevel[homologyType][elem1].append(elem2)
	else :
		print "Patate 6 : problem"
	return(DicoHomologyLevel[homologyType][elem1])
#----------------------------------------------------------------------#
def ImportHomologySp1(l, currentSp, DicoHomology, dicoTypeIdSpecies, specie1, specie2):
	""" We fill the DicoHomology with the information from the file of the specie 1
	Parameters
	    ----------
	     l : string, one line of the inputfile
	     currentSp : string, correspond to the current sp of the file
	     DicoHomology : {"Gene": {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								  "Paralogy Sp1" : {idGeneSp1 : [idGeneSp1]}
								  "Paralogy Sp2" : {idGeneSp2 : [idGeneSp2]}}
						"Transcript" : {"Orthology" : {idTranscririptSp1 : [idTranscriptSp2]}}
										"Paralogy Sp1" : {{idTranscririptSp1 : [idTranscriptSp1]}}
										"Paralogy Sp2" : {{idTranscririptSp2 : [idTranscriptSp2]}} 
	    Returns
	    -------
	     DicoHomology : like in the parameters but with new informations
	"""
	Specie1GeneID, Specie1TranscriptID, Specie2GeneID, Specie2TranscriptID = GetIDTypeSpecie(dicoTypeIdSpecies, specie1, specie2)
	
	GeneSpecie1, TranscriptSpecie1, GeneSpecie1Paralogue, TranscriptSpecie1Paralogue, GeneSpecie2, TranscriptSpecie2 = InitialisationHomology(l, currentSp)
	
	DicoHomology["Gene"]["Orthologues"][GeneSpecie1] = AddGeneHomology(GeneSpecie1, GeneSpecie2, DicoHomology["Gene"], "Orthologues")
	DicoHomology["Gene"]["Paralogues Sp1"][GeneSpecie1] = AddGeneHomology(GeneSpecie1, GeneSpecie1Paralogue, DicoHomology["Gene"], "Paralogues Sp1")

	# same but for transcript
	if GeneSpecie2 and TranscriptSpecie2 : # we verify if there is some orthologues
		if re.search(Specie2TranscriptID, TranscriptSpecie2): # first we need to verify
			DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] = AddGeneHomology(TranscriptSpecie1, TranscriptSpecie2, DicoHomology["Transcript"], "Orthologues")
			DicoHomology["Transcript"]["Paralogues Sp1"][TranscriptSpecie1] = AddGeneHomology(TranscriptSpecie1, TranscriptSpecie1Paralogue, DicoHomology["Transcript"], "Paralogues Sp1")
		else: # the id is the id of a protein
			if TranscriptSpecie2 in DicoGeneTranscriptProteinSpecie2["Protein"] :
				# if the prot id is in the dico GTP so get the id of the transcript
				TranscFromProtSpecie2 = DicoGeneTranscriptProteinSpecie2["Protein"][TranscriptSpecie2][0]
				DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] = AddGeneHomology(TranscriptSpecie1, TranscFromProtSpecie2, DicoHomology["Transcript"], "Orthologues")
	if GeneSpecie1Paralogue and TranscriptSpecie1Paralogue : # we verify if there is some paralogues
		if re.search(Specie1TranscriptID, TranscriptSpecie1Paralogue):
			#if the transcript for the paralogy is a protein, we change it
			if TranscriptSpecie1Paralogue in DicoGeneTranscriptProteinSpecie2["Protein"] :
				TranscFromProtSpecie1 = DicoGeneTranscriptProteinSpecie2["Protein"][TranscriptSpecie1Paralogue][0]
				DicoHomology["Transcript"]["Paralogues Sp1"][TranscriptSpecie1] = AddGeneHomology(TranscriptSpecie1, TranscFromProtSpecie1, DicoHomology["Transcript"], "Paralogues Sp1")
	return(DicoHomology)
#----------------------------------------------------------------------#
def ImportHomologySp2(l, currentSp, DicoHomology, dicoTypeIdSpecies, specie1, specie2, dicoGTP):
	""" When the DicoHomology is already created and fill with the information 
		of the specie 1, we check if there is some update with the file of the specie 2
	Parameters
	    ----------
	     l : string, one line of the inputfile
	     currentSp : string, correspond to the current sp of the file
	     DicoHomology : {"Gene": {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								  "Paralogy Sp1" : {idGeneSp1 : [idGeneSp1]}
								  "Paralogy Sp2" : {idGeneSp2 : [idGeneSp2]}}
						"Transcript" : {"Orthology" : {idTranscririptSp1 : [idTranscriptSp2]}}
										"Paralogy Sp1" : {{idTranscririptSp1 : [idTranscriptSp1]}}
										"Paralogy Sp2" : {{idTranscririptSp2 : [idTranscriptSp2]}} 
	    Returns
	    -------
	     DicoHomology : like in the parameters but with new informations
	"""
	Specie1GeneID, Specie1TranscriptID, Specie2GeneID, Specie2TranscriptID = GetIDTypeSpecie(dicoTypeIdSpecies, specie1, specie2)
	
	GeneSpecie1, TranscriptSpecie1, GeneSpecie2Paralogue, TranscriptSpecie2Paralogue, GeneSpecie2, TranscriptSpecie2 = InitialisationHomology(l, currentSp)
	
	DicoHomology["Gene"]["Paralogues Sp2"][GeneSpecie1] = AddGeneHomology(GeneSpecie1, GeneSpecie2Paralogue, DicoHomology["Gene"], "Paralogues Sp2")
	
	DicoHomology.update(MajOrthology(GeneSpecie1, TranscriptSpecie1, GeneSpecie2, TranscriptSpecie2, DicoHomology))
	
	# same but for transcript
	if GeneSpecie2 and TranscriptSpecie2 : # we verify if there is some orthologues
		if re.search(Specie2TranscriptID, TranscriptSpecie2): # first we need to verify
			DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] = AddGeneHomology(TranscriptSpecie1, TranscriptSpecie2, DicoHomology["Transcript"], "Orthologues")
			DicoHomology["Transcript"]["Paralogues Sp2"][TranscriptSpecie1] = AddGeneHomology(TranscriptSpecie1, TranscriptSpecie2Paralogue, DicoHomology["Transcript"], "Paralogues Sp2")
		else: # the id is the id of a protein
			if TranscriptSpecie2 in dicoGTP["Protein"] :
				# if the prot id is in the dico GTP so get the id of the transcript
				TranscFromProtSpecie2 = dicoGTP["Protein"][TranscriptSpecie2][0]
				DicoHomology["Transcript"]["Orthologues"][TranscriptSpecie1] = AddGeneHomology(TranscriptSpecie1, TranscFromProtSpecie2, DicoHomology["Transcript"], "Orthologues")
	if GeneSpecie2Paralogue and TranscriptSpecie2Paralogue : # we verify if there is some paralogues
		if not re.search(Specie1TranscriptID, TranscriptSpecie2Paralogue):
			#if the transcript for the paralogy is a protein, we change it
			if TranscriptSpecie2Paralogue in dicoGTP["Protein"] :
				TranscFromProtSpecie2 = dicoGTP["Protein"][TranscriptSpecie2Paralogue][0]
				DicoHomology["Transcript"]["Paralogues Sp2"][TranscriptSpecie2] = AddGeneHomology(TranscriptSpecie2, TranscFromProtSpecie2, DicoHomology["Transcript"], "Paralogues Sp2")
	return(DicoHomology)
#----------------------------------------------------------------------#
def ImportHomology(filename, DicoGeneTranscriptProteinSpecie2, dicoTypeIdSpecies, specie1, specie2, currentSp, DicoHomology):
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
	print "Import of Homology " + str(currentSp) + " : "
	
	cpt = 0
	
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			cpt += 1
			if l :
				if currentSp == "Sp1" :
					DicoHomology.update(ImportHomologySp1(l, currentSp, DicoHomology, dicoTypeIdSpecies, specie1, specie2))
				else:
					DicoHomology.update(ImportHomologySp2(l, currentSp, DicoHomology, dicoTypeIdSpecies, specie1, specie2, DicoGeneTranscriptProteinSpecie2))
				
			stat = round((float(cpt)/float(len(lines)))*100,1)
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
	return(DicoHomology)
#----------------------------------------------------------------------#
def importInfoG4Transcript(filenameG4InTranscript, DicoInfoG4):
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
					DicoInfoG4["Transcript"][G4][idTranscript]["Localisations"] = [localisation]
					DicoInfoG4["Transcript"][G4][idTranscript]["Biotype"] = [biotype]
				elif G4 in DicoInfoG4["Transcript"] and idTranscript not in DicoInfoG4["Transcript"][G4] :
					# the G4 is encounter but not in this transcript
					DicoInfoG4["Transcript"][G4][idTranscript] = {"Localisations" : []}
					DicoInfoG4["Transcript"][G4][idTranscript]["Localisations"] = [localisation]
					DicoInfoG4["Transcript"][G4][idTranscript]["Biotype"] = [biotype]
	return(DicoInfoG4)
#----------------------------------------------------------------------#
def importInfoG4Gene(DicoGTP, DicoInfoG4):
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
			if Transcrit in DicoGTP["Transcript"] : 
				gene = DicoGTP["Transcript"][Transcrit].keys() # get only the name of the transcript's gene
			else:
				gene = [""]
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
				if len(G4PerTranscript[Transcrit]["Localisations"]) > 1 or len(G4PerTranscript[Transcrit]["Biotype"]) > 1 :
					print "PATATE"
				localisationTranscript = G4PerTranscript[Transcrit]["Localisations"][0]
				biotypeTranscript = G4PerTranscript[Transcrit]["Biotype"][0]
				# 2 if because there can those two condition can be verified
				if localisationTranscript not in localisationGene :
					DicoInfoG4["Gene"][G4][gene]["Localisations"].append(localisationTranscript)
				if biotypeTranscript not in biotypeGene : 
					DicoInfoG4["Gene"][G4][gene]["Biotype"].append(biotypeTranscript)
	return(DicoInfoG4)
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
	
	DicoInfoG4.update(importInfoG4Transcript(filenameG4InTranscript, DicoInfoG4))
	DicoInfoG4.update(importInfoG4Gene(DicoGTP, DicoInfoG4))
	
	return(DicoInfoG4)
