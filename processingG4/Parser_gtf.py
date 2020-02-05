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

def addIntrons(df):
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
				cptIntron += 1
	dfIntron = dfIntron.reset_index(drop=True)
	dfIntron['Coords'] = [ [dfIntron.Start[x], dfIntron.End[x]] for x in range(0,len(dfIntron))]
	return(dfIntron)

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
	dfTmp = pd.DataFrame()
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('exon') ].dropna())
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('CDS') ].dropna())
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('five_prime_utr') ].dropna())
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('three_prime_utr') ].dropna())
	dfTmp = dfTmp.append(df[ df.Feature.str.contains('UTR') ].dropna())
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
	dfTmp = dfTmp.append( addIntrons(dfTmp) )
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

def importUTRFromGTF(filename):
	"""Imports gtf file.

	This function aims to retrieve as more as possible informations contained
	in the gtf file related to transcripts.

	:param filelame: name of the gtf file.
	:type filename: string

	:returns: dicoTr, contains all transcripts from a psecie and all features
		linked to it : UTR, exon, intron.
	:rtype: dictionary
	"""
	exists = os.path.isfile(filename)
	if exists :
		dicoTr = {}
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					words = l.split('\t')
					attributes = words[8].split(';')
					idGene = attributes[0].split('"')[1]
					feature = words[2]
					chrm = words[0]
					startFeature = int(words[3])
					endFeature = int(words[4])
					strand = words[6]
					strand = changeStrandFormat(strand)
					biotype = retrieveBiotypeFronAttributes(attributes, feature)
					if feature == 'five_prime_utr' :
						idTr = retrieveIdTrFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, '5UTR'))
						dicoTr[idTr]['5UTR'].update({'Chromosome' : chrm,
													'Start' : startFeature,
													'End' :endFeature,
													'Biotype' : biotype,
													'Strand' : strand})
					elif feature == 'three_prime_utr' :
						idTr = retrieveIdTrFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, '3UTR'))
						dicoTr[idTr]['3UTR'].update({'Chromosome' : chrm,
													'Start' : startFeature,
													'End' :endFeature,
													'Biotype' : biotype,
													'Strand' : strand})
				elif re.search('genome-version', l):
					assembly = l.split(' ')[1]
	else:
		print("This file don't exist : " + filename)
	return dicoTr

def importGTF(filename):
	"""Imports gtf file.

	This function aims to retrieve as more as possible informations contained
	in the gtf file related to transcripts.

	:param filelame: name of the gtf file.
	:type filename: string

	:returns: dicoTr, contains all transcripts from a psecie and all features
		linked to it : UTR, exon, intron.
	:rtype: dictionary
	"""
	exists = os.path.isfile(filename)
	if exists :
		dicoTr = {}
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					words = l.split('\t')
					attributes = words[8].split(';')
					idGene = attributes[0].split('"')[1]
					feature = words[2]
					chrm = words[0]
					startFeature = int(words[3])
					endFeature = int(words[4])
					strand = words[6]
					strand = changeStrandFormat(strand)
					biotype = retrieveBiotypeFronAttributes(attributes, feature)
					if feature == 'exon' :
						idTr = retrieveIdTrFronAttributes(attributes)
						rank, idExon = retrieveInfoExonFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, 'Exon'))
						dicoTr[idTr]['Exon'].update({idExon : {'Chromosome' : chrm,
																'Start' : startFeature,
																'End' : endFeature,
																'Biotype' : biotype,
																'Strand' : strand,
																'Rank' : rank}})
					elif feature == 'five_prime_utr' :
						idTr = retrieveIdTrFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, '5UTR'))
						dicoTr[idTr]['5UTR'].update({'Chromosome' : chrm,
													'Start' : startFeature,
													'End' :endFeature,
													'Biotype' : biotype,
													'Strand' : strand})
					elif feature == 'three_prime_utr' :
						idTr = retrieveIdTrFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, '3UTR'))
						dicoTr[idTr]['3UTR'].update({'Chromosome' : chrm,
													'Start' : startFeature,
													'End' :endFeature,
													'Biotype' : biotype,
													'Strand' : strand})
					elif feature == 'transcript':
						idTr = retrieveIdTrFronAttributes(attributes)
						if idTr not in dicoTr:
							dicoTr.update({idTr : {}})
						dicoTr[idTr].update({'Gene' : idGene,
											'Chromosome' : chrm,
											'Start' : startFeature,
											'End' : endFeature,
											'Biotype' : biotype,
											'Strand' : strand,
											'Assembly' : assembly})
				elif re.search('genome-version', l):
					assembly = l.split(' ')[1]
	else:
		print("This file don't exist : " + filename)
	return dicoTr

def importLocationFronGTF(filename):
	"""Imports gtf file.

	This function aims to retrieve informations of locations (Exon, CDS, UTR,
	start and end codon).

	:param filelame: name of the gtf file.
	:type filename: string

	:returns: dicoTr, contains all transcripts from a psecie and all features
		linked to it : UTR, exon, intron.
	:rtype: dictionary
	"""
	exists = os.path.isfile(filename)
	dico = {}
	if exists :
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			listLoc = ['exon', 'CDS', 'five_prime_utr', 'start_codon',
				'stop_codon', 'three_prime_utr']
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					words = l.split('\t')
					attributes = words[8].split(';')
					gene = attributes[0].split('"')[1]
					feature = words[2]
					chr = words[0]
					startFeature = int(words[3])
					endFeature = int(words[4])
					strand = words[6]
					strand = changeStrandFormat(strand)
					biotype = retrieveBiotypeFronAttributes(attributes, 'transcript')
					if feature in listLoc:
						idTr = retrieveIdTrFronAttributes(attributes)
						if feature == 'stop_codon' or feature == 'start_codon':
							startFeature = startFeature - 50
							endFeature = endFeature + 50
						idLoc = chr +':'+ str(startFeature) +'~'+ \
							str(endFeature) +':'+ strand
						if chr not in dico:
							dico[chr] = {gene : { feature : {idLoc : []} } }
						if gene not in dico[chr]:
							dico[chr][gene] = { feature : {idLoc : []} }
						if feature not in dico[chr][gene]:
							dico[chr][gene][feature] = {idLoc : []}
						if idLoc not in dico[chr][gene][feature] :
							dico[chr][gene][feature][idLoc] = []
						dico[chr][gene][feature][idLoc].append( idTr+'-'+biotype )
	else:
		print("This file don't exist : " + filename)
	return dico

def computeLength(filename):
	dicoGene = importGTFGene(filename)
	geneLength = 0
	dfTr = importGTFdf(filename)
	trLength = 0
	for gene in dicoGene:
		geneLength += dicoGene[gene]['End'] - dicoGene[gene]['Start'] +1
	transcripts = list(set(dfTr['Transcript']))
	for tr in transcripts:
		dfTmp = dfTr[ dfTr.Transcript == tr ].dropna().reset_index(drop=True)
		if (dfTmp['Type'][0] == 'Non Coding' and
			('five_prime_utr' in dfTmp['Feature'] or
			'three_prime_utr' in dfTmp['Feature'])):
			pass
		else:
			start = min(dfTmp['Start'])
			end = max(dfTmp['End'])
			trLength += end - start +1
	# dfTmp = dfTr[ dfTr.Feature == 'Exon' ].dropna()
	# starts = dfTmp['End']
	# counts = Counter(starts)
	# dupids = [id for id in starts if counts[id] > 1]
	# print dupids
	return geneLength, trLength

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp = arg.specie	# specie to parse
	print(sp)
	filename = "/home/anais/Documents/Data/Genomes/" + sp + \
		"/" + sp + ".gtf"
	df = importGTFdf(filename, sp)
	# print(len(list(set(df.Gene))))
	# print(len(list(set(df.Transcript))))
	# dfTmp = pd.DataFrame()
	# dfTmp = dfTmp.append(df[ df.Type == 'Coding'])
	# print(len(list(set(dfTmp.Transcript))))
	# dfTmp = pd.DataFrame()
	# dfTmp = dfTmp.append(df[ df.Type == 'Non coding'])
	# print(len(list(set(dfTmp.Transcript))))
	# dfTmp = pd.DataFrame()
	# dfTmp = dfTmp.append(df[ df.Feature == 'five_prime_utr'])
	# print(len(list(set(dfTmp.Transcript))))
	# dfTmp = pd.DataFrame()
	# dfTmp = dfTmp.append(df[ df.Feature == 'three_prime_utr'])
	# print(len(list(set(dfTmp.Transcript))))
	# gene, tr = computeLength(filename)
	# print('Gene length : ' + str(gene))
	# print('Transcript length : ' + str(tr))
