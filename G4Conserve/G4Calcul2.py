#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
Copyright Université de Sherbrooke, Département de Biochimie et Département d'Informatique

contact: sarah.belhamiti@usherbrooke.ca

This software is a computer program whose calcul the statistic value of PG4r (density and percent)

---------------------------------------------------

``G4Calcul`` **module description**:

From a tabulate files which contain G4 regions detected in one chromosome (or all chromosomes), this module will study the the distribution of the regions in the transcriptome. 

.. moduleauthor:: Sarah.Belhamiti

"""

import sys
import os
import re
import argparse
from pprint import pprint
import recurentFunction as rF

def GetLengthFraction(positionA, positionB):
	""" This fonction is to define the lenght between two positions 
	(positionA and positionB). It Doesn't matter the order of position 
	(positionA can be > positionB)
	    Parameters
	    ----------
	    positionA : string
	    positionB : string
	    Returns
	    -------
	    length: integer
	"""
	length = 0
	if ((positionA and positionB) != ""): # because some transcript doesn't have 5' for exemple (and start5 and end5 will be == '')
		length = (abs(int(positionA)-int(positionB)))+1
	return int(length)

def GetLengthCodant(exonList):
	""" This fonction is to define the lenght of exonique part.Can be use for the intronique part.
	    Parameters
	    ----------
	    exonList : list
		list of couple exon contain in one transcript
	    Returns
	    -------
	    lengthCodant: integer
		lenght of exonique part
	"""
	lengthCodant = 0
	if ('' not in exonList): # if transcript contain exon
		for coupleExon in exonList: # for each exon
			startExon = coupleExon.split('-')[0] 
			endExon = coupleExon.split('-')[1]
			lengthCodant += (GetLengthFraction(startExon,endExon))# add the lenght of this exon to the total (lengthCodant)
	return lengthCodant 

def GetLengthTotal(exonList):
	""" This fonction is to define the lenght total of a transcript based on the list of exon
	    Parameters
	    ----------
	    exonList : list
		list of couple exon contain in one transcript
	    Returns
	    -------
	    lengthCodant: integer
		lenght of exonique part
	"""
	lengthTotal = 0
	allExonsBorder = [] # all borders of exons contained in a transcript
	for coupleExon in exonList:
		allExonsBorder.append(int(coupleExon.split('-')[0]))
		allExonsBorder.append(int(coupleExon.split('-')[1]))
	lengthTotal=GetLengthFraction(min(allExonsBorder),max(allExonsBorder))
	return lengthTotal

def readLineIndex(line):
	words = line.split('|')
	dicoLine = {'idTr' : words[0],
				'idGene' : words[1],
				'Chromosome' : words[2],
				'Strand' : words[3],
				'geneBiotype' : words[4],
				'exonList' : words[5].split(';'),
				'intronList' : words[6].split(';'),
				'start5' : words[7],
				'end5' : words[8],
				'start3' : words[9],
				'end3' : words[10].rstrip()}
	return dicoLine


def GetLenghtByTranscript(filename, BiotypeByTranscript, AnnotationTranscript,Coding ):
	""" Create dictionary of lenght for each transcript
	if transcript is coding ,format =  [length 5,length cds,length 3,lenght exon,lenght intron,number junction 5_codant,number junction codant_Intron,number junction Intron_codant,number junction codant_3,length total, number junction CDS_CDS ]
	else format = [lenght exon,lenght intron,number junction codant_Intron,number junction Intron_codant,length total, number junction CDS_CDS ]
	    Parameters
	    ----------
	     filename : string
		name of the file which contain informatiob fr each transcript (except the transcript's biotype)
 	     BiotypeByTranscript : dico
		dictionary, where key =id of transcript and value = biotype of this transcript
	    AnnotationTranscript: dico
		dictionary with information of annotation for the transcripts By Ensembl (True if good anotation, else False)
	    Coding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    Returns
	    -------
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	"""
	LenghtByTranscript = {}
	inputfile = open(filename,"r") # file opening for reading
	for line in inputfile: # for each transcript
		dicoLine = readLineIndex(line)
		annotationTranscript = AnnotationTranscript.get(dicoLine["idTr"])
		transcriptBiotype = BiotypeByTranscript.get(dicoLine["idTr"])
		if annotationTranscript and transcriptBiotype :
			dicoTr = {"Biotype" : transcriptBiotype}
			if transcriptBiotype in Coding:
				dicoTr["Length 5UTR"] = GetLengthFraction(dicoLine["start5"],
														dicoLine["end5"])
				dicoTr["Length exons"] = GetLengthCodant(dicoLine["exonList"])
				dicoTr["Length 3UTR"] = GetLengthFraction(dicoLine["start3"],
														dicoLine["end3"])
				dicoTr["Length CDS"] = dicoTr["Length exons"] - \
										dicoTr["Length 5UTR"] - \
										dicoTr["Length 3UTR"]
				if dicoLine["intronList"][0] != '':
					dicoTr["Length introns"] = GetLengthCodant(dicoLine["intronList"])
				dicoTr["Number overlap 5UTR-coding"] = 1
				dicoTr["Number overlap coding-intron"] = len(dicoLine["exonList"]) - 1
				dicoTr["Number overlap intron-coding"] = len(dicoLine["exonList"]) - 1
				dicoTr["Number overlap coding-UTR3"] = 1
				dicoTr["Total length"] = GetLengthTotal(dicoLine["exonList"])
				dicoTr["Number junction"] = len(dicoLine["exonList"]) - 1
				dicoTr["Number of splicingSites"] = dicoTr["Number overlap coding-intron"] + \
													dicoTr["Number overlap intron-coding"] + \
													dicoTr["Number junction"]
				#(number junction on splicing events = CDS_Intron + INtronCDS+CDS_CDS)
			else:
				dicoTr["Length exons"] = GetLengthCodant(dicoLine["exonList"])
				if dicoTr["Length exons"] != '':
					dicoTr["Length introns"] = GetLengthCodant(dicoLine["intronList"])
				dicoTr["Number overlap coding-intron"] = len(dicoLine["exonList"]) - 1
				dicoTr["Number overlap intron-coding"] = len(dicoLine["exonList"]) - 1
				dicoTr["Total length"] = GetLengthTotal(dicoLine["exonList"])
				dicoTr["Number junction"] = len(dicoLine["exonList"]) - 1
				dicoTr["Number of splicingSites"] = dicoTr["Number overlap coding-intron"] + \
													dicoTr["Number overlap intron-coding"] + \
													dicoTr["Number junction"]
			LenghtByTranscript[ dicoLine["idTr"] ] = dicoTr
	inputfile.close()
	return LenghtByTranscript

def iniDicoNumberG4Coding(biotype):
	dico = {"5" : 0,
			"CDS" : 0,
			"3" : 0,
			"Exon" : 0,
			"Intron" : 0,
			"junction_5_CDS": 0,
			"junction_CDS_Intron" : 0,
			"junction_Intron_CDS" : 0,
			"junction_CDS_3" : 0,
			"Total" : 0,
			"junction_CDS_CDS" : 0,
			"Splicing" : 0,
			"Biotype" : biotype}
	return dico

def iniDicoNumberG4NonCoding(biotype):
	dico = {"ExonNC" : 0,
			"IntronNC" : 0,
			"junction_ExonNC_IntronNC" : 0,
			"junction_IntronNC_ExonNC" : 0,
			"junction_ExonNC_ExonNC" : 0,
			"Total" : 0,
			"Splicing" : 0,
			"Biotype" : biotype}
	return dico

def readLineG4InTr(line):
	words = line.split('\t')
	dicoLine = {"InfoG4ByTr" : words[0],
				"idTr" : words[0].split('|')[0],
				"Chromosome" : words[0].split('|')[1].split(':')[0],
				"g4Start" : int(words[0].split('|')[1].split(':')[1].split('-')[0]),
				"g4End" : int(words[0].split('|')[1].split(':')[1].split('-')[1]),
				"Strand" : words[0].split('|')[2],
				"Location" : words[5],
				"Biotype" : words[6].rstrip()}
	return dicoLine

def GetNumbersG4Transcript(filename, Coding):
	""" Create dictionary of number of G4 for each transcript depending the localisation
	if transcript is coding ,format =  [number G4 in section 5, cds, 3,exon, intron,junction 5, junction CDS, junction i, junction 3, total, junction cds_cds ]
	else format = [number G4 in exon, intron, junction CDS, junction i, total, junction cds_cds ]
	    Parameters
	    ----------
	     filename : string
		name of the file which contain information for each G4 in a transcript 
	    Coding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    Returns
	    -------
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	"""
	NumbersG4Transcript={}
	inputfile= open(filename,"r") # file opening for reading
	for line in inputfile: # for each transcript
		if ('InfoG4ByTranscript' not in line and line != '\n'):
			dicoL = readLineG4InTr(line)
			if dicoL["Biotype"] in Coding:
				if NumbersG4Transcript.has_key(dicoL["idTr"]) == False:
					numberG4 = iniDicoNumberG4Coding(dicoL["Biotype"])
				else: # if gene contain G4 detected 
					numberG4 = NumbersG4Transcript.get(dicoL["idTr"])
				numberG4[ dicoL["Location"] ] += 1
				numberG4["Total"] += 1
				if (dicoL["Location"] == '5' or dicoL["Location"] == 'CDS' or 
					dicoL["Location"] == '3' or dicoL["Location"] == 'junction_5_CDS'
					or dicoL["Location"] == 'junction_CDS_3'):
					numberG4["Exon"] += 1
				if (dicoL["Location"] == 'junction_CDS_CDS' or
					dicoL["Location"] == 'junction_Intron_CDS' or
					dicoL["Location"] == 'junction_CDS_Intron'):
					numberG4["Splicing"] += 1
				NumbersG4Transcript[ dicoL["idTr"] ] = numberG4
			else:
				if (NumbersG4Transcript.has_key(dicoL["idTr"]) == False):
					numberG4 = iniDicoNumberG4NonCoding(dicoL["Biotype"])
				else:
					numberG4 = NumbersG4Transcript.get(dicoL["idTr"])
				numberG4[ dicoL["Location"] ] += 1
				NumbersG4Transcript[ dicoL["idTr"] ] = numberG4
				numberG4["Total"] += 1
				if (dicoL["Location"] == 'junction_ExonNC_ExonNC' or 
					dicoL["Location"] == 'junction_IntronNC_ExonNC' or
					dicoL["Location"] == 'junction_ExonNC_IntronNC'):
					numberG4["Splicing"] += 1
				NumbersG4Transcript[dicoL["idTr"]] = numberG4
	inputfile.close()
	return NumbersG4Transcript

def iniListLoca(biotype, coding):
	if biotype in coding:
		localisationName = ['5', 'CDS', '3', 'Exon', 'Intron',
							'junction_5_CDS', 'junction_CDS_3',
							'junction_CDS_Intron', 'junction_Intron_CDS',
							'Total', 'junction_CDS_CDS','Splicing']
	else:
		localisationName = ['ExonNC', 'IntronNC', 'junction_ExonNC_IntronNC',
							'junction_IntronNC_ExonNC', 'Total', 'junction_ExonNC_ExonNC',
							'Splicing']
	return localisationName

def GetNumberTranscriptByCat(dicoLength, NumbersG4Transcript, coding):
	""" Create dictionary of number of transcript in each categorie 
	    Parameters
	    ----------
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 
	    Coding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    Returns
	    -------
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	"""
	NumberTrByCat = {}
	for idTr in dicoLength:
		pprint(dicoLength[idTr])
		if (NumberTrByCat.has_key(dicoLength[idTr]["Biotype"]) == False):
			# if biotype not contain in the dico
			localisationName = iniListLoca(dicoLength[idTr]["Biotype"], coding)
			if len(localisationName) == 12: # coding
				dicoTr = iniDicoNumberG4Coding(dicoLength[idTr]["Biotype"])
			elif len(localisationName) == 7: # non coding
				dicoTranscript = iniDicoNumberG4Coding(dicoLength[idTr]["Biotype"])
			dicoTr["# tr in biotype"] = 0
		else: # if biotype already contained in the dico
			dicoTr = NumberTrByCat.get(dicoLength[idTr]["Biotype"])
		dicoTr["# tr in biotype"] += 1
		if NumbersG4Transcript.has_key(idTr) == True:
			numberG4 = NumbersG4Transcript.get(idTr)
			pprint(numberG4)
			for position in numberG4:
				if position != "Biotype" and numberG4[position] != 0:
					# if the transcript contain G4 at this position
					dicoTr[position] += dicoTr[position]
		NumberTrByCat[ dicoLength[idTr]["Biotype"] ] = dicoTr
	return NumberTrByCat

def CreateRnaTypeDetected(LenghtByTranscript, rnaTypeDetected):
	""" Create list of RNA type where PG4r were detected 
	    Parameters
	    ----------
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     rnaTypeDetected : list
		list of RNA type where PG4r were detected
	    Returns
	    -------
	     rnaTypeDetected : list
		list of RNA type where PG4r were detected
	"""
	rnaTypeDetected=[]
	for key, values in LenghtByTranscript.items():
		biotypeTranscript=key.split('-')[1]
		if (biotypeTranscript not in rnaTypeDetected) :
			rnaTypeDetected.append(biotypeTranscript)
	return rnaTypeDetected

def EnrichssmentByFamily(familyName, family, state, 
						LenghtByTranscript, NumbersG4Transcript,
						NumberTranscriptByCat):
	""" Create list of calcul for a complete RNA family
	    Parameters
	    ----------
	    familyName : string
		Family where PG4r was detectd
	    family : list
		list of type contained in this family
	    state : string
		condirion where PG4r is detected (~ location)
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 

	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	if (familyName== 'Coding'):
		localisationName=['5','CDS','3','Exon','Intron','junction_5_CDS','junction_CDS_3','junction_CDS_Intron','junction_Intron_CDS','Total','junction_CDS_CDS']
	else:
		localisationName=['ExonNC','IntronNC','junction_ExonNC_IntronNC','junction_IntronNC_ExonNC','Total','junction_ExonNC_ExonNC']
	position=GetPositionForLocalisation(state, familyName)
	length=0
	num=0
	listeBiotype=[]
	for key, values in LenghtByTranscript.items(): # for each transcript in the transcriptom
		transcriptId=key.split('-')[0]
		transcriptBiotype=key.split('-')[1]
		if (transcriptBiotype in family): # if transcript in this family
			length=length+values[position]
			if (transcriptBiotype not in listeBiotype):
				listeBiotype.append(transcriptBiotype)
			if (NumbersG4Transcript.has_key(key)==True): # if transcript contain G4
				num=num+NumbersG4Transcript.get(key)[position]
	numTrAll=0
	numTrWithG4=0
	for biotype in listeBiotype:
		position=GetPositionForLocalisation(state, familyName)
		numTrAll=numTrAll+NumberTranscriptByCat.get(biotype)[0]
		numTrWithG4=numTrWithG4+NumberTranscriptByCat.get(biotype)[position+1]
	percent=round(numTrWithG4/float(numTrAll)*100,2)
	if ('junction' in state):
		return [familyName,length, num, round(num/float(length),4),numTrAll,numTrWithG4,percent]
	else:
		return [familyName,length, num, round(num/float(length)*1000,4),numTrAll,numTrWithG4,percent]

def EnrichssmentByRnaType(familyName, family, state,
						LenghtByTranscript, NumbersG4Transcript,
						rnaTypeDetected, NumberTranscriptByCat):
	""" Create list of calcul for a RNA type  
	    Parameters
	    ----------
	    familyName : string
		Family where PG4r was detectd
	    family : list
		list of type contained in this family
	    state : string
		condirion where PG4r is detected (~ location)
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	     rnaTypeDetected : list
		list of RNA type where PG4r were detected
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 

	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	if (familyName== 'Coding'):
		localisationName=['5','CDS','3','Exon','Intron','junction_5_CDS','junction_CDS_3','junction_CDS_Intron','junction_Intron_CDS','Total','junction_CDS_CDS']
	else:
		localisationName=['ExonNC','IntronNC','junction_ExonNC_IntronNC','junction_IntronNC_ExonNC','Total','junction_ExonNC_ExonNC']
	position=GetPositionForLocalisation(state, familyName)
	liste=[]
	for rnaType in family:
		if (rnaType in rnaTypeDetected):
			length=0
			num=0
			numTrAll=0
			numTrWithG4=0
			for key, values in LenghtByTranscript.items(): # for each transcript in the transcriptom
				transcriptId=key.split('-')[0]
				transcriptBiotype=key.split('-')[1]
				if (transcriptBiotype == rnaType): # if transcript in this family
					length=length+values[position]
					if (NumbersG4Transcript.has_key(key)==True): # if transcript contain G4
						num=num+NumbersG4Transcript.get(key)[position]
			position=GetPositionForLocalisation(state,familyName)
			numTrAll=NumberTranscriptByCat.get(rnaType)[0]
			numTrWithG4=NumberTranscriptByCat.get(rnaType)[position+1]
			if (length !=0):
				percent=round(numTrWithG4/float(numTrAll)*100,2)
				if ('junction' in state):
					liste.append([rnaType,length, num, round(num/float(length),4),numTrAll,numTrWithG4,percent])
				else:
					liste.append([rnaType,length, num, round(num/float(length)*1000,4),numTrAll,numTrWithG4,percent])
			else:
					liste.append([rnaType,length, num, round(num/float(length)*1000,4),numTrAll,numTrWithG4,0])
			
	return liste

def EnrichssmentByRnaTypeLocalisation (familyName, family,state,LenghtByTranscript,NumbersG4Transcript, rnaType,NumberTranscriptByCat):	
	""" Create list of calcul for a RNA type  at a specific position 
	    Parameters
	    ----------
	    familyName : string
		Family where PG4r was detectd
	    family : list
		list of type contained in this family
	    state : string
		condirion where PG4r is detected (~ location)
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	    rnaTypeDetected : list
		list of RNA type where PG4r were detected
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 
	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	position=GetPositionForLocalisation(state, familyName)
	length=0
	num=0
	numTrAll=0
	numTrWithG4=0
	for key, values in LenghtByTranscript.items(): # for each transcript in the transcriptom
		transcriptId=key.split('-')[0]
		transcriptBiotype=key.split('-')[1]
		if (transcriptBiotype == rnaType): # if transcript in this family
			length=length+values[position]
			if (NumbersG4Transcript.has_key(key)==True): # if transcript contain G4
				num=num+NumbersG4Transcript.get(key)[position]

	position=GetPositionForLocalisation(state,familyName)
	numTrAll=NumberTranscriptByCat.get(rnaType)[0]
	numTrWithG4=NumberTranscriptByCat.get(rnaType)[position+1]

	percent=round(numTrWithG4/float(numTrAll)*100,2)

	if (length !=0):
		if ('junction' in state):
			return [rnaType,length, num, round(num/float(length),4),numTrAll,numTrWithG4,percent]
		else:
			return [rnaType,length, num, round(num/float(length)*1000,4),numTrAll,numTrWithG4,percent]
	else:
		return [rnaType,length, num, 0,numTrAll,numTrWithG4,percent]

def EnrichssmentByFamilyLocalisation(familyName, family,state,LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat):
	""" Create ensembl of list of calcul for all RNA type  at a specific position for a category
	    Parameters
	    ----------
	    familyName : string
		Family where PG4r was detectd
	    family : list
		list of type contained in this family
	    state : string
		condirion where PG4r is detected (~ location)
	    LenghtByTranscript : dictionary
		dictionary of lenght for each transcript
	     NumbersG4Transcript : dictionary
		dictionary of number of G4 for each transcript depending the localisation
	    rnaTypeDetected : list
		list of RNA type where PG4r were detected
	    NumberTranscriptByCat :  dictionary
		dictionary of number of transcript for each category 
	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	liste=[]
	position=GetPositionForLocalisation(state, familyName)
	for rnaType in family:
		if (rnaType in rnaTypeDetected):
			liste.append(EnrichssmentByRnaTypeLocalisation (familyName, family,state,LenghtByTranscript,NumbersG4Transcript, rnaType,NumberTranscriptByCat ))
	if ('junction' in state):
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste])),4)])
	else:
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste])),4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)]) # *100 because percent
	return liste

def PrintInFile (outfilename,header,column,liste):
	""" print line for each element in a liste
	    Parameters
	    ----------
	    outfilename : string
		name file output
	    header : string
		header of the table
	    column : list
		list of column name
	    liste : list
		list of informations
	    Returns
	    -------
	    liste :  list
		list of calcul
	"""
	outfilename.write(header+'\n')
	outfilename.write(column+'\n')
	for element in liste:
		s=[str(item) for item in element]
		outfilename.write('\t'.join(s)+'\n')
	outfilename.write('\n')

def createDicoFamily():
	dicoFam = {"Coding" : ['IG_C_gene', 'IG_D_gene', 'IG_J_gene',
							'IG_LV_gene', 'IG_M_gene', 'IG_V_gene',
							'IG_Z_gene', 'nonsense_mediated_decay',
							'nontranslating_CDS', 'non_stop_decay',
							'protein_coding', 'TR_C_gene', 'TR_D_gene',
							'TR_gene', 'TR_J_gene', 'TR_V_gene'],
				"Pseudogene" : ['transcribed_unitary_pseudogene',
								'disrupted_domain', 'IG_C_pseudogene',
								'IG_J_pseudogene', 'IG_pseudogene',
								'IG_V_pseudogene', 'processed_pseudogene',
								'pseudogene', 'transcribed_processed_pseudogene',
								'transcribed_unprocessed_pseudogene',
								'translated_processed_pseudogene',
								'translated_unprocessed_pseudogene',
								'TR_J_pseudogene', 'TR_V_pseudogene',
								'unitary_pseudogene', 'unprocessed_pseudogene',
								'polymorphic_pseudogene'],
				"Long non codint" : ['macro_lncRNA', 'bidirectional_promoter_lncRNA',
									'sense_intronic', '3prime_overlapping_ncRNA',
									'ambiguous_orf', 'antisense',
									'lincRNA', 'ncrna_host','non_coding',
									'processed_transcript', 'retained_intron',
									'sense_overlapping'],
				"Short non coding" : ['vaultRNA', 'scaRNA', 'miRNA',
										'miRNA_pseudogene', 'misc_RNA', 
										'misc_RNA_pseudogene', 'Mt_rRNA',
										'Mt_tRNA', 'Mt_tRNA_pseudogene',
										'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA',
										'RNase_P_RNA', 'rRNA', 'rRNA_pseudogene', 
										'scRNA', 'scRNA_pseudogene', 'snlRNA',
										'snoRNA', 'snoRNA_pseudogene', 'snRNA',
										'snRNA_pseudogene', 'SRP_RNA', 'tmRNA',
										'tRNA', 'tRNA_pseudogene','ribozyme'],
				"Predicted" : ['TEC'],
				'Non coding' : ['vaultRNA', 'macro_lncRNA',
								'bidirectional_promoter_lncRNA',
								'sense_intronic', '3prime_overlapping_ncRNA',
								'ambiguous_orf', 'antisense', 'lincRNA',
								'ncrna_host', 'non_coding', 
								'processed_transcript', 'retained_intron',
								'sense_overlapping', '3prime_overlapping_ncRNA',
								'scaRNA', 'miRNA', 'miRNA_pseudogene', 'misc_RNA',
								'misc_RNA_pseudogene', 'Mt_rRNA', 'Mt_tRNA',
								'Mt_tRNA_pseudogene', 'ncRNA', 'pre_miRNA',
								'RNase_MRP_RNA', 'RNase_P_RNA', 'rRNA',
								'rRNA_pseudogene', 'scRNA', 'scRNA_pseudogene',
								'snlRNA', 'snoRNA', 'snoRNA_pseudogene',
								'snRNA', 'snRNA_pseudogene', 'SRP_RNA',
								'tmRNA', 'tRNA', 'tRNA_pseudogene', 'ribozyme']}
	return dicoFam

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	GITDIR=os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR+'results/')
	parser.add_argument ('-chr', '--chromosome', default = 'all')
	parser.add_argument ('-specie', '--specie', default = 'HS')
	parser.add_argument ('-c', '--choice', default = 'Coding')
	return parser

def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path # directory which contain all the directory chromosome
	chromosome = arg.chromosome
	specie = arg.specie
	choice = arg.choice
	GITDIR=os.getcwd()+'/'

	Coding=['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene', 'IG_M_gene', 'IG_V_gene', 'IG_Z_gene', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay', 'protein_coding', 'TR_C_gene', 'TR_D_gene', 'TR_gene', 'TR_J_gene', 'TR_V_gene']
	Pseudogene=['transcribed_unitary_pseudogene','disrupted_domain', 'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene', 'IG_V_pseudogene', 'processed_pseudogene', 'pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene', 'translated_unprocessed_pseudogene', 'TR_J_pseudogene', 'TR_V_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene','polymorphic_pseudogene']
	LongNC=['macro_lncRNA','bidirectional_promoter_lncRNA','sense_intronic','3prime_overlapping_ncRNA','ambiguous_orf','antisense','lincRNA','ncrna_host','non_coding','processed_transcript','retained_intron','sense_overlapping']
	ShortNC=['vaultRNA','scaRNA','miRNA', 'miRNA_pseudogene', 'misc_RNA', 'misc_RNA_pseudogene', 'Mt_rRNA' ,'Mt_tRNA', 'Mt_tRNA_pseudogene', 'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA', 'RNase_P_RNA', 'rRNA', 'rRNA_pseudogene','scRNA', 'scRNA_pseudogene', 'snlRNA', 'snoRNA', 'snoRNA_pseudogene', 'snRNA', 'snRNA_pseudogene', 'SRP_RNA', 'tmRNA', 'tRNA', 'tRNA_pseudogene','ribozyme']
	Predictif=['TEC']
	NonCoding=['vaultRNA','macro_lncRNA','bidirectional_promoter_lncRNA','sense_intronic','3prime_overlapping_ncRNA','ambiguous_orf','antisense','lincRNA','ncrna_host','non_coding','processed_transcript','retained_intron','sense_overlapping','3prime_overlapping_ncRNA','scaRNA','miRNA', 'miRNA_pseudogene', 'misc_RNA', 'misc_RNA_pseudogene', 'Mt_rRNA' ,'Mt_tRNA', 'Mt_tRNA_pseudogene', 'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA', 'RNase_P_RNA', 'rRNA', 'rRNA_pseudogene','scRNA','scRNA_pseudogene', 'snlRNA', 'snoRNA', 'snoRNA_pseudogene', 'snRNA', 'snRNA_pseudogene', 'SRP_RNA', 'tmRNA', 'tRNA', 'tRNA_pseudogene','ribozyme']
	SupraFamily=[Coding, NonCoding, Pseudogene, Predictif]
	SupraFamilyName=['Coding', 'NonCoding', 'Pseudogene', 'Predictif']
	dicoFam = createDicoFamily()
	
	if (chromosome == 'all'):
		directory = path+'/all'
		index = directory+'/'+specie+'_transcript_unspliced_All_Index.txt'	
		indextranscriptBiotype = directory+'/transcriptType_All'
		fileG4InTranscriptome = directory+'/HS_All_G4InTranscript.txt'
	else:
		directory = path+'/chr'+chromosome
		index = directory+'/'+specie+'_transcript_unspliced_chr'+chromosome+'_Index.txt'
		indextranscriptBiotype = path+'/transcriptType/transcriptType_chr'+chromosome
		fileG4InTranscriptome = '/home/local/USHERBROOKE/bels2814/Documents/These/HUMAN/HS_chr'+chromosome+'_G4InTranscript.txt'

	# Assignation of variable for the start
	NumberTranscriptByCat={} # dictionary of list of number of each for each rna type (position 0 = all gene, position 1 = gene with G4)
	rnaTypeDetected=[] # list of rnaType detected in this species
	
	# Filling dictionary
	BiotypeByTranscript = rF.createDictionaryBiotypeByTranscript(indextranscriptBiotype)
	AnnotationTranscript = rF.GetAnnotationTranscript(index, dicoFam["Coding"],
													BiotypeByTranscript)
	LenghtByTranscript = GetLenghtByTranscript(index,
												BiotypeByTranscript,
												AnnotationTranscript,
												dicoFam["Coding"])
	NumbersG4Transcript = GetNumbersG4Transcript(fileG4InTranscriptome,
												dicoFam["Coding"])
	NumberTranscriptByCat = GetNumberTranscriptByCat(LenghtByTranscript,
													NumbersG4Transcript,
													dicoFam["Coding"])
	
	rnaTypeDetected = CreateRnaTypeDetected(LenghtByTranscript,rnaTypeDetected)
	
	column='RNA Type\tTotal Length\t#PG4r\tDensity\t#Tr\t#TrWithPG4r\t%Tr'
	outfilename= open(GITDIR+'resultsTable/'+choice+'.csv',"w")
	
	if (choice == 'Coding'):
		liste=[]
		condition='Total'
		liste.append(EnrichssmentByFamily ('Coding', Coding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('NonCoding', NonCoding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Pseudogene', Pseudogene,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Predictif', Predictif,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header='Whole transcriptome'
		PrintInFile (outfilename,header,column,liste)
		###
		liste=EnrichssmentByRnaType ('Coding', Coding,'Total',LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header='By mRNA Subclass '
		PrintInFile (outfilename,header,column,liste)
		###
		conditions=[['Exon','Intron','Splicing'],['5','CDS','3', 'junction_5_CDS', 'junction_CDS_3'],['junction_CDS_Intron','junction_Intron_CDS','junction_CDS_CDS']]
		headers=[['Exon','Intron','Splicing Site'],["5'UTR",'CDS',"3'UTR", 'Codon start', 'Codon stop'],['Donor splicing site','Acceptor splicing site','Splicing junction']]
		for condition in conditions:
			for state in condition:
				#liste=EnrichssmentByFamily ('Coding', Coding,state,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat)
				liste=EnrichssmentByFamilyLocalisation('Coding', Coding,state,LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
				header=headers[conditions.index(condition)][condition.index(state)]
				PrintInFile (outfilename,header,column,liste)

	elif (choice == 'NonCoding'):
		liste=[]
		condition='Total'
		liste.append(EnrichssmentByFamily ('Coding', Coding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('NonCoding', NonCoding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Pseudogene', Pseudogene,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Predictif', Predictif,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header='Whole transcriptome'
		PrintInFile (outfilename,header,column,liste)
		###
		liste.append(EnrichssmentByFamily ('LongNC', LongNC,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('ShortNC', ShortNC,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header='By ncRNA groupe'
		PrintInFile (outfilename,header,column,liste)
		###
		liste=EnrichssmentByRnaType ('LongNC', LongNC,'Total',LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header=header='By long ncRNA Subclass '
		PrintInFile (outfilename,header,column,liste)
		###
		conditions=[['ExonNC','IntronNC','Splicing'],['junction_ExonNC_IntronNC','junction_IntronNC_ExonNC','junction_ExonNC_ExonNC']]
		headers=[['Exon','Intron','Splicing Site'],['Donor splice site','Acceptor splice site','Splicing Junction ']]
		for condition in conditions:
			for state in condition:
				#liste=EnrichssmentByFamily ('LongNC', LongNC,state,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat)
				liste=EnrichssmentByFamilyLocalisation('LongNC', LongNC,state,LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
				header=headers[conditions.index(condition)][condition.index(state)]
				PrintInFile (outfilename,header,column,liste)
		###
		liste=EnrichssmentByRnaType ('ShortNC', ShortNC,'Total',LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])
		header=header='By short ncRNA Subclass '
		PrintInFile (outfilename,header,column,liste)
	
	
	else:
		

		liste=[]
		condition='Total'
		liste.append(EnrichssmentByFamily ('Coding', Coding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('NonCoding', NonCoding,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Pseudogene', Pseudogene,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(EnrichssmentByFamily ('Predictif', Predictif,condition,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat))
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])

		header='Whole transcriptome'
		PrintInFile (outfilename,header,column,liste)
		
		
		#print '-----------------------------------------------------'


		liste=EnrichssmentByRnaType ('Pseudogene', Pseudogene,'Total',LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
		liste.append(['Overall',sum([x[1]for x in liste]),sum([x[2]for x in liste]),round(sum([x[2]for x in liste])/float(sum([x[1]for x in liste]))*1000,4),sum([x[4]for x in liste]),sum([x[5]for x in liste]),round(sum([x[5]for x in liste])/float(sum([x[4]for x in liste]))*100,4)])



		header='By pseudogene Subclass '
		PrintInFile (outfilename,header,column,liste)
		####
		conditions=[['ExonNC','IntronNC','Splicing'],['junction_ExonNC_IntronNC','junction_IntronNC_ExonNC','junction_ExonNC_ExonNC']]
		headers=[['Exon','Intron','Splicing Site'],['Donor splice site','Acceptor splice site','Splicing Junction ']]
		for condition in conditions:
			
			for state in condition:
				#liste=EnrichssmentByFamily ('Pseudogene', Pseudogene,state,LenghtByTranscript,NumbersG4Transcript,NumberTranscriptByCat)
				
				liste=EnrichssmentByFamilyLocalisation('Pseudogene', Pseudogene,state,LenghtByTranscript,NumbersG4Transcript, rnaTypeDetected,NumberTranscriptByCat)
				header=headers[conditions.index(condition)][condition.index(state)]
				PrintInFile (outfilename,header,column,liste)
	


main()
