#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	January 2020

Description:
	This script will read ouputs from G4RNA screener for the shuffled dataset.
		This script is different from the WT beacause sequences are annotated
		differently.

Command Line:
	For one chr : python ~/PATH/getMainDensities.py
"""

import os
import re
import argparse
import pandas as pd
from pprint import pprint
import recurrentFunction as rF

def getInfo(df):
	"""Retrieves informations of a windows and parse it into a dictionary.

	As gene windows and junction windows are not formated the same way, this
	function aims to parse them into the same type of dictionary.

	:param df: contain all overlaping windows.
	:type df: dataFrame

	:returns: dico, contains all infromation for one window.
	:rtype: dictionary
	"""
	geneDesc = df.geneDesc.iloc[0]
	lastRow = len(df.index) - 1
	dico = {'geneId' : geneDesc,
			'lastRow' : len(df.index) - 1,
			'meancGcC' : df.cGcC.mean(),
			'meanG4H' : df.G4H.mean(),
			'meanG4NN' : df.G4NN.mean(),
			'pG4Start' : int(df.wStart.iloc[0]),
			'pG4End' : int(df.wEnd.iloc[lastRow])}
	return dico

def mergeOverlappingSequences(dfTmp):
	"""Merge the sequences of overlaping windows.

	:param dfTmp: contain overlaping windows.
	:type dfTmp: dataFrame

	:returns: seq, sequence merged.
	:rtype: string
	"""
	seq = str(dfTmp.seqG4.iloc[0])
	dfTmp = dfTmp.sort_values(by=['wStart'])
	for w in range(1,len(dfTmp)):
		step = int(dfTmp.wStart.iloc[w] - dfTmp.wStart.iloc[w-1])
		# convert to int elsewise it's a float
		wSeq = dfTmp.seqG4.iloc[w]
		seq += wSeq[-step:]
	return seq

def mergeWindows(df, junctionLength):
	"""Merge overlaping windows.

	:param df: contain overlaping windows.
	:type df: dataFrame
	:param junctionLength: length of unction, by default it's 100 nt.
	:type junctionLength: integer

	:returns: pG4, contains the pG4 which is the merge of overlaping windows.
	:rtype: dictionary
	"""
	pG4rSeq = mergeOverlappingSequences(df)
	dicoInfo = getInfo(df)
	pG4Start = dicoInfo['pG4Start']
	pG4End = dicoInfo['pG4End']
	pG4 = {'id' : [ dicoInfo['geneId'] ],
			'cGcC' : [ dicoInfo['meancGcC'] ],
			'G4H' : [ dicoInfo['meanG4H'] ],
			'G4NN' : [ dicoInfo['meanG4NN'] ],
			'Start' : [pG4Start], 'End' : [pG4End],
			'seqG4' : [pG4rSeq]}
	return pG4

def mergeG4(df, dicoParam):
	"""Browses all junction window to find those that are overlapping.

	Here we browse all junctions windows. We will only kept those that overlap
	the 100 nucleotid. Indeed, if the window over thresholds don't overlap this
	position, it only in a gene and not a junction.

	:param df: contains all windows.
	:type df: dataFrame
	:param dicoParam: contains all parameters that were given to g4rna screener.
	:type dicoParam: dictionnary

	:returns: dfpG4, contain all pG4 for that strand.
	:rtype: dataFrame
	"""
	dfTmp = pd.DataFrame()
	dfpG4 = pd.DataFrame()
	dfTmp = dfTmp.append(df[0:1]) # store the first window
	if len(df) == 1:
		dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
				dicoParam["junctionLength"]))
		dfpG4 = dfpG4.append(dfTmp)
	else:
		for index, row in df.iterrows(): # w for window
			# browses all windows over thresholds, exept the first one
			if (row.geneDesc == df.geneDesc.iloc[index-1] and
				((row.wStart >= df.wStart.iloc[index-1] and row.wStart <= df.wEnd.iloc[index-1]) or
				(row.wEnd >= df.wStart.iloc[index-1] and row.wEnd <= df.wEnd.iloc[index-1]))):
				# if window overlap, add window at the current pG4
				dfTmp = dfTmp.append(df[index:index+1])
				if index == len(df)-1 :
					dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
							dicoParam["junctionLength"]))
					dfpG4 = dfpG4.append(dfTmp)
			else: # new pG4
				dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
						dicoParam["junctionLength"]))
				dfpG4 = dfpG4.append(dfTmp)
				dfTmp = df.iloc[index:index+1]
				if index == len(df)-1 :
					dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
							dicoParam["junctionLength"]))
					dfpG4 = dfpG4.append(dfTmp)
				# reinitiate the dfTmp with the new pG4
	return dfpG4

def filterOnScores(dicoParam, dfWindows):
	"""Filter the windows based on thresholds.

	:param dicoParam: contains all parameters that were given to g4rna screener.
	:type dicoParam: dictionnary
	:param dfWindows: contains all windows of all genes from one specie.
	:type dfWindows: dataframe

	:returns: dfWindows, with only windows upper thresholds.
	:rtype: dataFrame
	"""
	dfWindows = dfWindows[ dfWindows.cGcC >= dicoParam["cGcC"] ].dropna()
	dfWindows = dfWindows[ dfWindows.G4H >= dicoParam["G4H"] ].dropna()
	dfWindows = dfWindows[ dfWindows.G4NN >= dicoParam["G4NN"] ].dropna()
	return dfWindows

def getpG4(filename, dicoParam):
	dfpG4 = pd.DataFrame()
	try:
		dfWindows = pd.read_csv(filename, sep='\t', index_col=0)
	except:
		print("This file couldn't be converted in data frame : " + filename)
	else:
		# dataFrame with all windows from G4RNA Screener
		dfWindows.columns = ['geneDesc','cGcC',
							'G4H','seqG4','wStart',
							'wEnd', 'G4NN']
		dfWindows = filterOnScores(dicoParam, dfWindows)
		dfWindows = dfWindows.reset_index(drop=True)
		dfpG4 = dfpG4.append(mergeG4(dfWindows, dicoParam))
		return dfpG4

def annotationpG4(dfpG4unNN, dicoType):
	dfpG4Ann = pd.DataFrame()
	dfpG4unNN.reset_index(drop=True, inplace=True)
	dicotmp = dfpG4unNN.to_dict('index')
	for row in dicotmp:
		listTrBt = dicotmp[row]['id'].split(':')[5].split('|')
		location = [dicotmp[row]['id'].split(':')[1]]
		coords = dicotmp[row]['id'].split(':')[2:5]
		for TrBt in listTrBt:
			tr = TrBt.split('-')[0]
			bt = TrBt.split('-')[1]
			if bt in dicoType:
				type = dicoType[bt]
			else:
				type = 'None'
			tmpRow = {'Transcript' : tr,
				'Location' : location, 'Sequence' : dicotmp[row]['seqG4'],
				'Start' : coords[1].split('~')[0],
				'End' : coords[1].split('~')[1], 'Strand' : coords[2],
				'cGcC' : dicotmp[row]['cGcC'], 'G4H' : dicotmp[row]['G4H'],
				'G4NN' : dicotmp[row]['G4NN'], 'Biotype' : bt, 'Type' : type}
			dfTmp = pd.DataFrame.from_dict(tmpRow)
			dfpG4Ann = dfpG4Ann.append(dfTmp)
	return dfpG4Ann

def mergeWindow(path, dicoParam, sp):
	directory = path + sp + '/ShuffleCSV'
	dfpG4 = pd.DataFrame()
	dicoType = rF.createDicoType()
	# directory containing data for a specie
	for path, dirs, files in os.walk(directory):
		# for each element of the directory to passed
		for filename in files: # for each files
			inputfile = directory + '/' + filename
			if ('location' in filename):
				# windows in genes
				unAnnopG4 = getpG4(inputfile, dicoParam)
				annopG4 = annotationpG4(unAnnopG4, dicoType)
				dfpG4 = dfpG4.append(annopG4)
				dfpG4 = dfpG4.reset_index(drop=True)
	if len(dfpG4) > 0:
		dfpG4 = dfpG4.drop_duplicates(subset=None, keep='first', inplace=False)
		dfpG4 = dfpG4.reset_index(drop=True)
	return dfpG4

def removeTr(path, df, loca):
	filename = path + 'listTrBm.txt'
	filtereddf = pd.DataFrame()
	with open(filename) as f:
		lines = f.read().splitlines()
		for index, row in df.iterrows():
			if loca == 'Intron':
				if len(row.id.split(':')) > 1:
					tr = row.id.split(':')[1]
				else:
					tr = ''
			else:
				tr = row.id.split(':')[2]
			if tr in lines:
				filtereddf = filtereddf.append(row)
	return filtereddf

def main(dicoParam, path, sp):
	dfpG4 = mergeWindow(path, dicoParam, sp)
	print ('\t'+str(dfpG4.shape))
	outputFN = path + sp + '/pG4_shuffled.csv'
	dfpG4.to_csv(path_or_buf = outputFN, header=True, index=None, sep='\t')

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
	parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-E', '--junctionLength', default = 100)
	parser.add_argument ('-W', '--WINDOW', default = 60)
	parser.add_argument ('-S', '--STEP', default = 10)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	sp = arg.specie
	dicoParam = rF.createDicoParam(arg)
	main(dicoParam, path, sp)
