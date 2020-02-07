#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	May 2019

Description:
	This script is principaly a library containing many methods to import a gtf
	file. It can be imported as a dictionary or as a dataFrame. While imported
	as dictionary, intron are computed and added to the dictionnary. An ouput
	containing introns information can be created.

Data availability:
	* gtf file via Ensembl ftp.

"""

import re
import os
import argparse
import numpy as np
import pandas as pd
from pprint import pprint
from collections import Counter
import recurrentFunction as rF

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_gtf')
	parser.add_argument ('-sp', '--specie', default = 'mus_musculus')
	return parser

def changeStrandFormat(strand):
	"""Changes the format of the strand from +/- to 1/-1.
	"""
	if strand == '+':
		strand = '1'
	elif strand == '-':
		strand = '-1'
	return strand

def retrieveBiotypeFronAttributes(attributes, feature):
	"""Gets the biotype from attributes.

	:param feature: transcript, gene, exon, or utr.
	:type feature: string
	:param attributes: last colons of gtf file, contains a lot of informations
		but is different depending on the feature.
	:type attributes: list

	:returns: biotype of the feature. The biotype of a gene can be different
		compared to the transcript biotype.
	:rtype: string
	"""
	biotype = ''
	for attribute in attributes:
		if re.search(feature+"_biotype", attribute):
			biotype = attribute.split('"')[1]
	return biotype

def retrieveIdTrFronAttributes(attributes):
	"""Gets the id transcript from attributes.

	:param attributes: last colons of gtf file, contains a lot of informations
		but is different depending on the feature.
	:type attributes: list

	:returns: idTr, transctipt identifier.
	:rtype: string
	"""
	idTr = ''
	for attribute in attributes:
		if re.search('transcript_id', attribute):
			idTr = attribute.split('"')[1]
	return idTr

def retrieveInfoExonFronAttributes(attributes):
	"""Gets rank and exon number from attributes.

	:param attributes: last colons of gtf file, contains a lot of informations
		but is different depending on the feature.
	:type attributes: list

	:returns: rank, idExon of an exon.
	:rtype: strings
	"""
	rank = 0
	idExon = 0
	for attribute in attributes:
		if re.search('exon_number', attribute):
			rank = int(attribute.split('"')[1])
		elif re.search('exon_id', attribute):
			idExon = attribute.split('"')[1]
	return rank, idExon

def retrieveUTRFromDico(dico):
	"""Gets UTR coords if they exist.

	:param dico: contains all transcripts and their features.
	:type dico: dictionary

	:returns: start5UTR, end5UTR, start3UTR, end3UTR.
	:rtypes: strings
	"""
	start5UTR = ''
	end5UTR = ''
	start3UTR = ''
	end3UTR = ''
	if '5UTR' in dico :
		start5UTR, end5UTR = retrieveCoordUTRFromDico(dico['5UTR'])
	if '3UTR' in dico :
		start3UTR, end3UTR = retrieveCoordUTRFromDico(dico['3UTR'])
	return start5UTR, end5UTR, start3UTR, end3UTR

def retrieveCoordUTRFromDico(dico):
	start = dico["Start"]
	end = dico["End"]
	return start, end

def createKeyTranscript(dico, idGene, idTr, feature):
	if idTr not in dico :
		dico[idTr] = {"Gene" : idGene,
						feature : {}}
	if feature not in dico[idTr] :
		dico[idTr].update({feature : {}})
	return dico

def importTranscriptFromGTF(filename):
	dicoTr = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: # browse all lines
			if not l.startswith('#') and l:
				words=l.split('\t')
				attributes = words[8].split(';')
				idGene = attributes[0].split('"')[1]
				feature = words[2]
				chrm = words[0]
				startFeature = words[3]
				endFeature = words[4]
				strand = words[6]
				strand = changeStrandFormat(strand)
				biotype = retrieveBiotypeFronAttributes(attributes, feature)
				if feature == "exon" :
					idTr = retrieveIdTrFronAttributes(attributes)
					rank, idExon = retrieveInfoExonFronAttributes(attributes)
					dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, "Exon"))
					dicoTr[idTr]["Exon"].update({idExon : {"Chromosome" : chrm,
															"Start" : int(startFeature),
															"End" :int(endFeature),
															"Biotype" : biotype,
															"Strand" : strand,
															"Rank" : rank}})
				elif feature == "transcript":
					idTr = retrieveIdTrFronAttributes(attributes)
					if idTr not in dicoTr:
						dicoTr.update({idTr : {}})
					dicoTr[idTr].update({"Gene" : idGene,
										"Chromosome" : chrm,
										"Start" : int(startFeature),
										"End" : int(endFeature),
										"Biotype" : biotype,
										"Strand" : strand,
										"Type" : addTypeTr(biotype)})
	return dicoTr

def importGTFGene(filename):
	"""Imports genes features from the gtf file.

	:param filename: name of the gtf file.
	:type filename: string

	:returns: dicoGene, contains all informations abotu genes.
	:rtype: dictionary
	"""
	dicoGene = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: # browse all lines
			if not l.startswith('#') and l:
				words = l.split('\t')
				geneId = words[8].split(';')[0].split('"')[1]
				attributes = words[8].split(';')
				feature = words[2]
				if feature == "gene" :
					chrm = words[0]
					startFeature = int(words[3])
					endFeature = int(words[4])
					strand = words[6]
					biotype = retrieveBiotypeFronAttributes(attributes, feature)
					if strand == '+':
						strand = '1'
					elif strand == '-':
						strand = '-1'
					if geneId not in dicoGene:
						dicoGene[geneId] = {}
					dicoGene[geneId].update({'Chromosome' : chrm,
											'Start' : startFeature,
											'End' : endFeature,
											'Strand' : strand,
											'Biotype' : biotype,
											'Assembly' : assembly})
				if feature == 'transcript':
					geneId = attributes[0].split('"')[1]
					idTr = retrieveIdTrFronAttributes(attributes)
					if geneId not in dicoGene:
						dicoGene[geneId] = {'Transcript' : []}
					elif 'Transcript' not in dicoGene[geneId]:
						dicoGene[geneId].update({'Transcript' : []})
					dicoGene[geneId]['Transcript'].append(idTr)
			elif re.search('genome-version', l):
				assembly = l.split(' ')[1]
	return dicoGene

def extractDicoGeneTr(filename):
	dicoGene = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: # browse all lines
			if not l.startswith('#') and l:
				words = l.split('\t')
				geneId = words[8].split(';')[0].split('"')[1]
				attributes = words[8].split(';')
				feature = words[2]
				if feature == 'transcript':
					geneId = attributes[0].split('"')[1]
					idTr = retrieveIdTrFronAttributes(attributes)
					if geneId not in dicoGene:
						dicoGene[geneId] = []
					dicoGene[geneId].append(idTr)
	return(dicoGene)

def getDicoTrGene(filename):
	dicoGene = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: # browse all lines
			if not l.startswith('#') and l:
				words = l.split('\t')
				attributes = words[8].split(';')
				feature = words[2]
				if feature == 'transcript':
					geneId = attributes[0].split('"')[1]
					idTr = retrieveIdTrFronAttributes(attributes)
					biotype = retrieveBiotypeFronAttributes(attributes, feature)
					dicoGene[idTr] = {'Gene' : geneId, 'Biotype' : biotype}
	return(dicoGene)

def addTranscript(attributes):
	"""Apply function to retrieve the transcript id of a feature.

	:param attributes: all attributes from the gtf file, corresponds to the
		ninth column of the file.
	:type attributes: string

	:returns: transcript id.
	:rtype: string
	"""
	attributes = attributes.split(';')
	return retrieveIdTrFronAttributes(attributes)

def addGene(attributes):
	"""Apply function to retrieve the transcript id of a feature.

	:param attributes: all attributes from the gtf file, corresponds to the
		ninth column of the file.
	:type attributes: string

	:returns: transcript id.
	:rtype: string
	"""
	attributes = attributes.split(';')
	geneID = attributes[0].split('"')[1]
	return geneID

def addBiotype(attributes):
	"""Apply function to retrieve the biotype of a transcript.

	:param attributes: all attributes from the gtf file, corresponds to the
		ninth column of the file.
	:type attributes: string

	:returns: biotype of the transcript.
	:rtype: string
	"""
	attributes = attributes.split(';')
	return retrieveBiotypeFronAttributes(attributes, 'transcript')

def addTypeTr(biotype):
	"""Apply function to retrieve the type of a transcript depending on biotype.

	:param biotype: biotype of a transcript.
	:type biotype: string

	:returns: type, coding or non coding.
	:rtype: string
	"""
	coding = ['IG_C_gene', 'IG_D_gene', 'IG_J_gene',
			'IG_LV_gene', 'IG_M_gene', 'IG_V_gene',
			'IG_Z_gene', 'nonsense_mediated_decay', 'non_stop_decay',
			'protein_coding', 'TR_C_gene', 'TR_D_gene',
			'TR_gene', 'TR_J_gene', 'TR_V_gene']
	if biotype in coding:
		return 'Coding'
	else:
		return 'Non coding'

def mapUTRonCDS(coordUTR, coordCDS, strand):
	if ( (coordUTR[1] < coordCDS[0] and strand == '+') or
		(coordUTR[1] > coordCDS[0] and strand == '-') ):
		return 'five_prime_utr'
	else:
		return 'three_prime_utr'

def getUTR(df):
	dicoUTR = {}
	dfUTR = pd.DataFrame()
	dfCDS = pd.DataFrame()
	dfUTR = dfUTR.append( df[ df.Feature == 'UTR' ] )
	dfCDS = dfCDS.append( df[ df.Feature == 'CDS'])
	dicoTmp = dfUTR[['index1', 'Coords', 'Transcript', 'Strand']]
	dicoTmp = dicoTmp.set_index('index1').to_dict()
	dicoCDS = dfCDS[['Transcript', 'Coords']]
	dicoCDS = dicoCDS.set_index('Transcript').to_dict()
	for index in dicoTmp['Transcript']:
		CDScoords = pd.DataFrame()
		tr = dicoTmp['Transcript'][index]
		utr = mapUTRonCDS(dicoTmp['Coords'][index], dicoCDS['Coords'][tr], \
			dicoTmp['Strand'][index])
		dicoUTR[index] = utr
	df = df.replace({'index1' : dicoUTR})
	return df

def addIntrons(df, dicoTr):
	dfIntron = pd.DataFrame()
	dfExon = pd.DataFrame()
	dfExon = df[ df.Feature == 'exon']
	groups = dfExon.groupby('Transcript')
	for name, group in groups:
		if len(group) > 1:
			# if more then 1 exon, then there is an intron
			group = group.sort_values(by=['Start'])
			group = group.reset_index(drop=True)
			starts = group.Start
			ends = group.End
			cptIntron = 0
			# print(group.Coords)
			while cptIntron < len(starts) - 1:
				start = ends[cptIntron] + 1
				end = starts[cptIntron + 1] - 1
				if start != end and start < end:
					row = {'Feature' : ['intron'], 'Start' : start,
					'End' : end, 'Strand' : group.Strand[0],
					'Attributes' : group.Attributes[0],
					'Chromosome' : group.Chromosome[0],
					'Transcript' : group.Transcript[0], 'Gene' : group.Gene[0],
					'Biotype' : group.Biotype[0], 'Type' : group.Type[0]}
					dfRow = pd.DataFrame.from_dict(row)
					dfIntron = dfIntron.append(dfRow)
					if (start - 40 < dicoTr[ group.Transcript[0] ]['Start']):
						startN = dicoTr[ group.Transcript[0] ]['Start']
					else:
						startN -= 40
					if (end + 40 > dicoTr[ group.Transcript[0] ]['End']):
						endN = dicoTr[ group.Transcript[0] ]['End']
					else:
						endN += 40
					row = {'Feature' : ['junction'], 'Start' : startN,
					'End' : endN, 'Strand' : group.Strand[0],
					'Attributes' : group.Attributes[0],
					'Chromosome' : group.Chromosome[0],
					'Transcript' : group.Transcript[0], 'Gene' : group.Gene[0],
					'Biotype' : group.Biotype[0], 'Type' : group.Type[0]}
					dfRow = pd.DataFrame.from_dict(row)
					dfIntron = dfIntron.append(dfRow)
					if (start - 40 < dicoTr[ group.Transcript[0] ]['Start']):
						startN = dicoTr[ group.Transcript[0] ]['Start']
					else:
						startN -= 40
					if (start + 40 > dicoTr[ group.Transcript[0] ]['End']):
						endN = dicoTr[ group.Transcript[0] ]['End']
					else:
						endN += 40
					row = {'Feature' : ['donor'], 'Start' : startN,
					'End' : endN, 'Strand' : group.Strand[0],
					'Attributes' : group.Attributes[0],
					'Chromosome' : group.Chromosome[0],
					'Transcript' : group.Transcript[0], 'Gene' : group.Gene[0],
					'Biotype' : group.Biotype[0], 'Type' : group.Type[0]}
					dfRow = pd.DataFrame.from_dict(row)
					dfIntron = dfIntron.append(dfRow)
					if (end - 40 < dicoTr[ group.Transcript[0] ]['Start']):
						startN = dicoTr[ group.Transcript[0] ]['Start']
					else:
						startN -= 40
					if (end + 40 > dicoTr[ group.Transcript[0] ]['End']):
						endN = dicoTr[ group.Transcript[0] ]['End']
					else:
						endN += 40
					row = {'Feature' : ['acceptor'], 'Start' : startN,
					'End' : end, 'Strand' : group.Strand[0],
					'Attributes' : group.Attributes[0],
					'Chromosome' : group.Chromosome[0],
					'Transcript' : group.Transcript[0], 'Gene' : group.Gene[0],
					'Biotype' : group.Biotype[0], 'Type' : group.Type[0]}
					dfRow = pd.DataFrame.from_dict(row)
					dfIntron = dfIntron.append(dfRow)
				cptIntron += 1
	dfIntron = dfIntron.reset_index(drop=True)
	dfIntron['Coords'] = [ [dfIntron.Start[x], dfIntron.End[x]] for x in range(0,len(dfIntron))]
	return(dfIntron)

def correctCoords(df, dfTmp):
	dfTmpTr = pd.DataFrame()
	dfTmpLoc = pd.DataFrame()
	locations = ['junction', 'acceptor', 'donor', 'start_codon', 'stop_codon']
	for loc in locations:
		dfTmpLoc = dfTmpLoc.append(dfTmp[ dfTmp.Feature == loc ])
	dfTmpTr = dfTmpTr.append(df[ df.Feature.str.contains('transcript') ].dropna())
	dfTmpTr['Transcript'] = dfTmpTr.Attributes.apply(addTranscript)
	dicoTr = dfTmpTr.set_index('Transcript').to_dict('index')
	for index, row in dfTmpCodon.iterrows():
		if (row.Start - 40 < dicoTr[row.Transcript]['Start']):
			dfTmp['Start'].iloc[index] = dicoTr[row.Transcript]['Start']
		if (row.End + 40 > dicoTr[row.Transcript]['End']):
			dfTmp['End'].iloc[index] = dicoTr[row.Transcript]['End']
	return dfTmp

def addCodon(df, dicoTr):
	dfTmpCodon = pd.DataFrame()
	dfTmpCodon = dfTmpCodon.append(df[ df.Feature.str.contains('start_codon') ].dropna())
	dfTmpCodon = dfTmpCodon.append(df[ df.Feature.str.contains('stop_codon') ].dropna())
	dfTmpCodon['Transcript'] = dfTmpCodon.Attributes.apply(addTranscript)
	dfTmpCodon['Start'] = dfTmpCodon['Start'] - 40
	dfTmpCodon['End'] = dfTmpCodon['End'] + 40
	dfTmpCodon = dfTmpCodon.reset_index(drop=True)
	for index, row in dfTmpCodon.iterrows():
		if (row.Start - 40 < dicoTr[row.Transcript]['Start']):
			dfTmpCodon['Start'].iloc[index] = dicoTr[row.Transcript]['Start']
		if (row.End + 40 > dicoTr[row.Transcript]['End']):
			dfTmpCodon['End'].iloc[index] = dicoTr[row.Transcript]['End']
	return dfTmpCodon

def parseDF(df):
	"""Parses a dataframe to retrieve only usefull data for me.

	We only kept soem features : exon, CDS, 5UTR and 3UTR. From attributes we
	only kept : Transcript id, Biotype of the transcript and the type of the
	transcript (coding or not).

	:param df: contains all transcripts feature for a species but non parsed.
	:type df: dataFrame

	:returns: dfTmp, parsed df.
	:rtype: dataFrame
	"""
	dfTmpTr = pd.DataFrame()
	dfTmpTr = dfTmpTr.append(df[ df.Feature.str.contains('transcript') ].dropna())
	dfTmpTr['Transcript'] = dfTmpTr.Attributes.apply(addTranscript)
	dicoTr = dfTmpTr.set_index('Transcript').to_dict('index')
	dfTmp = pd.DataFrame()
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('exon') ].dropna())
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('CDS') ].dropna())
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('five_prime_utr') ].dropna())
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('three_prime_utr') ].dropna())
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('UTR') ].dropna())
	dfTmpCodon = addCodon(df, dicoTr)
	dfTmp = dfTmp.append(dfTmpCodon.dropna())
	dfTmp = dfTmp.reset_index(drop=True)
	dfTmp['Coords'] = [ [dfTmp.Start[x], dfTmp.End[x]] for x in range(0,len(dfTmp))]
	dfTmp['Transcript'] = dfTmp.Attributes.apply(addTranscript)
	dfTmp['Gene'] = dfTmp.Attributes.apply(addGene)
	dfTmp['index1'] = dfTmp.index
	if 'Biotype' not in df.columns.values:
		dfTmp['Biotype'] = dfTmp.Attributes.apply(addBiotype)
	else:
		dfTmp2 = pd.DataFrame()
		dfTmp2 = dfTmp2.append(df)
		dfTmp2['Transcript'] = dfTmp2.Attributes.apply(addTranscript)
		dfTmp2 = dfTmp2[['Transcript', 'Biotype']]
		dicoTmp = dfTmp2.set_index('Transcript').to_dict()
		dfTmp['Biotype'] = dfTmp['Transcript'].map(dicoTmp['Biotype'])
	if ('five_prime_utrf' not in set(df.Feature) and
		'three_prime_utr' not in set(df.Feature)):
		dfTmp = getUTR(dfTmp)
	dfTmp['Type'] = dfTmp.Biotype.apply(addTypeTr)
	dfTmp = dfTmp.append( addIntrons(dfTmp, dicoTr) )
	# dfTmp = correctCoords(df, dfTmp)
	return dfTmp

def importGTFdf(filename, sp):
	"""Imports a gtf file into a dataframe.

	Read the gtf file into a csv, then change the column names. NI goes for
	'Not Important', it's column I will not use so I delete them. I also
	retrieve only what I need from attributes and then dilate it.

	:param filename: name of the gtf file.
	:type filename: string

	:returns: df, contains all transcripts feature for a species.
	:rtype: dataFrame
	"""
	try:
		df = pd.read_csv(filename, sep='\t', index_col=0, skiprows=5)
	except:
		print("This file couldn't be converted in data frame : " + filename)
	else:
		# dataFrame with all windows from G4RNA Screener
		df.columns = ['Biotype', 'Feature','Start','End',
					'NI1', 'Strand', 'NI2', 'Attributes']
		del df['NI1']
		del df['NI2']
		if sp not in ['homo_sapiens', 'pan_troglodytes', 'mus_musculus', 'gallus_gallus']:
			del df['Biotype']
		df['Chromosome'] = df.index
		df = parseDF(df)
		del df['Attributes']
		return df

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp = arg.specie	# specie to parse
	print(sp)
	filename = "/home/anais/Documents/Data/Genomes/" + sp + \
		"/" + sp + ".gtf"
	df = importGTFdf(filename, sp)
