#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	February 2020

Description:
	This script imports a gtf file, and computes all locations used in this
	study (intron, codon start and stop, junction). Once those locations are
	computed, the dataFrame is exported as a csv to get a faster access for
	further scripts.
	
"""

import re
import os
import argparse
import numpy as np
import pandas as pd
from pprint import pprint

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

def addAcctepor(df,  dicoTr, start, end):
	if (end - 40 < dicoTr[ df.Transcript[0] ]['Start']):
		startN = dicoTr[ df.Transcript[0] ]['Start']
	else:
		startN = end -40
	if (end + 40 > dicoTr[ df.Transcript[0] ]['End']):
		endN = dicoTr[ df.Transcript[0] ]['End']
	else:
		endN = end + 40
	if df.Strand[0] == '+':
		row = {'Feature' : ['acceptor'], 'Start' : startN,
		'End' : end, 'Strand' : df.Strand[0],
		'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
		'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
		'Biotype' : df.Biotype[0], 'Type' : df.Type[0]}
	else:
		row = {'Feature' : ['donor'], 'Start' : startN,
		'End' : endN, 'Strand' : df.Strand[0],
		'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
		'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
		'Biotype' : df.Biotype[0], 'Type' : df.Type[0]}
	return row

def addDonor(df,  dicoTr, start, end):
	if (start - 40 < dicoTr[ df.Transcript[0] ]['Start']):
		startN = dicoTr[ df.Transcript[0] ]['Start']
	else:
		startN = start - 40
	if (start + 40 > dicoTr[ df.Transcript[0] ]['End']):
		endN = dicoTr[ df.Transcript[0] ]['End']
	else:
		endN = start + 40
	if df.Strand[0] == '+':
		row = {'Feature' : ['donor'], 'Start' : startN,
		'End' : endN, 'Strand' : df.Strand[0],
		'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
		'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
		'Biotype' : df.Biotype[0], 'Type' : df.Type[0]}
	else:
		row = {'Feature' : ['acceptor'], 'Start' : startN,
		'End' : endN, 'Strand' : df.Strand[0],
		'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
		'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
		'Biotype' : df.Biotype[0], 'Type' : df.Type[0]}
	return row

def addJunction(df, dicoTr, start, end):
	if (start - 40 < dicoTr[ df.Transcript[0] ]['Start']):
		startN = dicoTr[ df.Transcript[0] ]['Start']
	else:
		startN = start
	if (end + 40 > dicoTr[ df.Transcript[0] ]['End']):
		endN = dicoTr[ df.Transcript[0] ]['End']
	else:
		endN = end
	row = {'Feature' : ['junction'], 'Start' : str(startN)+'|'+str(start),
	'End' : str(end)+'|'+str(endN), 'Strand' : df.Strand[0],
	'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
	'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
	'Biotype' : df.Biotype[0], 'Type' : df.Type[0]}
	return row

def addIntron(df, start, end):
	row = {'Feature' : ['intron'], 'Start' : start,
	'End' : end, 'Strand' : df.Strand[0],
	'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
	'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
	'Biotype' : df.Biotype[0], 'Type' : df.Type[0]}
	return row

def getPointLocation(df, dicoTr):
	dfIntron = pd.DataFrame()
	dfExon = pd.DataFrame()
	dfExon = dfExon.append(df[ df.Feature == 'exon'])
	groups = dfExon.groupby('Transcript')
	for name, group in groups:
		if len(group) > 1:
			# if more then 1 exon, then there is an intron
			group = group.sort_values(by=['Start'])
			group = group.reset_index(drop=True)
			starts = group.Start
			ends = group.End
			cptIntron = 0
			while cptIntron < len(starts) - 1:
				start = ends[cptIntron] + 1
				end = starts[cptIntron + 1] - 1
				if start != end and start < end:
					row = addIntron(group, start, end)
					dfRow = pd.DataFrame.from_dict(row)
					dfIntron = dfIntron.append(dfRow)
					row = addJunction(group,  dicoTr, start, end)
					dfRow = pd.DataFrame.from_dict(row)
					dfIntron = dfIntron.append(dfRow)
					row = addDonor(group,  dicoTr, start, end)
					dfRow = pd.DataFrame.from_dict(row)
					dfIntron = dfIntron.append(dfRow)
					row = addAcctepor(group,  dicoTr, start, end)
					dfRow = pd.DataFrame.from_dict(row)
					dfIntron = dfIntron.append(dfRow)
				cptIntron += 1
	dfIntron = dfIntron.reset_index(drop=True)
	dfIntron['Coords'] = [ [dfIntron.Start[x], dfIntron.End[x]] for x in range(0,len(dfIntron))]
	return(dfIntron)

def addCodon(df, dicoTr):
	dfTmpCodon = pd.DataFrame()
	dfTmpCodon = dfTmpCodon.append(df[ df.Feature.str.contains('start_codon') ].dropna())
	dfTmpCodon = dfTmpCodon.append(df[ df.Feature.str.contains('stop_codon') ].dropna())
	dfTmpCodon['Transcript'] = dfTmpCodon.Attributes.apply(addTranscript)
	dfTmpCodon = dfTmpCodon.reset_index(drop=True)
	for index, row in dfTmpCodon.iterrows():
		if (row.Start - 40 < dicoTr[row.Transcript]['Start']):
			dfTmpCodon['Start'].iloc[index] = dicoTr[row.Transcript]['Start']
		else:
			dfTmpCodon['Start'].iloc[index] = row.Start - 40
		if (row.End + 40 > dicoTr[row.Transcript]['End']):
			dfTmpCodon['End'].iloc[index] = dicoTr[row.Transcript]['End']
		else:
			dfTmpCodon['End'].iloc[index] = row.End + 40
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
	del dfTmpTr
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
	dfTmp = dfTmp.append( getPointLocation(dfTmp, dicoTr) )
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

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_gtf')
	parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp = arg.specie	# specie to parse
	filename = "/home/anais/Documents/Data/Genomes/" + sp + \
		"/" + sp + ".gtf"
	output = "/home/anais/Documents/Data/Genomes/" + sp + \
		"/" + sp + ".csv"
	df = importGTFdf(filename, sp)
	df.to_csv(path_or_buf=output, header=True, index=None, sep='\t')
