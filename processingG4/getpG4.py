#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
``init.py`` **module description**:
This module has as input the output of G4RNAScreener and a dictionnary
with all the parameters that were used for G4RNA Screener.
This script filter the output of G4RNA Screener to only keep windows
over thresholds. Then all those windows are merged if they are
overlapping. Overlapping windows are pG4.
.. moduleauthor:: Anaìs Vannutelli, Michel Garant, Sarah Bellamiti and Aida Ouangraoua
March 2019
Université de Sherbrooke Canada
Laboratoty CoBiUS and Jean-Pierre Perreault
"""

import pandas as pd
import numpy as np
from pprint import pprint

def mergeOverlappingSequences(dfTmp):
	"""Merge the sequences of overlaping windows.

	:param dfTmp: contain overlaping windows.
	:type dfTmp: dataFrame

	:returns: seq, sequence merged.
	:rtype: string
	"""
	seq = str(dfTmp.seqG4.iloc[0])
	for w in range(1,len(dfTmp)):
		step = int(dfTmp.wStart.iloc[w] - dfTmp.wStart.iloc[w-1])
		# convert to int elsewise it's a float
		wSeq = dfTmp.seqG4.iloc[w]
		seq += wSeq[-step:]
	return seq

def getInfo(df, feature, option):
	"""Retrieves informations of a windows and parse it into a dictionary.

	As gene windows and junction windows are not formated the same way, this
	function aims to parse them into the same type of dictionary.

	:param df: contain all overlaping windows.
	:type df: dataFrame
	:param feature: junction or gene.
	:type feature: string

	:returns: dico, contains all infromation for one window.
	:rtype: dictionary
	"""
	geneDesc = df.geneDesc.iloc[0]
	lastRow = len(df.index) - 1
	if option == 'Annotation':
		if feature == 'Gene':
			geneDescSplit = geneDesc.split(' ')[2].split(':')
			dico = {'geneId' : geneDesc.split(' ')[0]}
		else:
			geneDescSplit = geneDesc.split(':')
			id = geneDescSplit[3] +':'+ geneDescSplit[4]
			dico = {'geneId' : id}
		dico.update({'Strand' : geneDescSplit[5],
					'Chromosome' : geneDescSplit[2],
					'geneStart' : int(geneDescSplit[3]), #for junction it's the
					'geneEnd' : int(geneDescSplit[4]), #intron start and end.
					'lastRow' : len(df.index) - 1,
					'geneDesc' : geneDesc,
					'meancGcC' : df.cGcC.mean(),
					'meanG4H' : df.G4H.mean(),
					'meanG4NN' : df.G4NN.mean(),
					'pG4Start' : int(df.wStart.iloc[0]),
					'pG4End' : int(df.wEnd.iloc[lastRow])})
	elif option == 'Random':
		dico = {'geneId' : geneDesc,
				'lastRow' : len(df.index) - 1,
				'meancGcC' : df.cGcC.mean(),
				'meanG4H' : df.G4H.mean(),
				'meanG4NN' : df.G4NN.mean(),
				'pG4Start' : int(df.wStart.iloc[0]),
				'pG4End' : int(df.wEnd.iloc[lastRow])}
	return dico

def mergeWindows(df, feature, junctionLength, option):
	"""Merge overlaping windows.

	:param df: contain overlaping windows.
	:type df: dataFrame
	:param feature: junction or gene.
	:type feature: string
	:param junctionLength: length of unction, by default it's 100 nt.
	:type junctionLength: integer

	:returns: pG4, contains the pG4 which is the merge of overlaping windows.
	:rtype: dictionary
	"""
	pG4rSeq = mergeOverlappingSequences(df)
	dicoInfo = getInfo(df, feature, option)
	if feature == "Junction":
		if (dicoInfo['pG4Start'] < junctionLength and
			dicoInfo['pG4End'] > junctionLength):
			pG4Start, pG4End = dicoInfo['pG4Start'], dicoInfo['pG4End']
			pG4 = True
		else:
			pG4 = None
	else:
		pG4Start = dicoInfo['pG4Start']
		pG4End = dicoInfo['pG4End']
		pG4 = True
	if pG4:
		if option == 'Annotation':
			pG4 = {'id' : [ dicoInfo['geneId'] ],
					'Strand' : [ dicoInfo['Strand'] ],
					'Chromosome' : [ dicoInfo['Chromosome'] ],
					'cGcC' : [ dicoInfo['meancGcC'] ],
					'G4H' : [ dicoInfo['meanG4H'] ],
					'G4NN' : [ dicoInfo['meanG4NN'] ],
					'Start' : [pG4Start], 'End' : [pG4End],
					'seqG4' : [pG4rSeq],
					'Feature' : [feature]}
		elif option == 'Random':
			pG4 = {'id' : [ dicoInfo['geneId'] ],
					'cGcC' : [ dicoInfo['meancGcC'] ],
					'G4H' : [ dicoInfo['meanG4H'] ],
					'G4NN' : [ dicoInfo['meanG4NN'] ],
					'Start' : [pG4Start], 'End' : [pG4End],
					'seqG4' : [pG4rSeq]}
	return pG4

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
	dfWindows = dfWindows[ dfWindows.G4H >= dicoParam["g4H"] ].dropna()
	dfWindows = dfWindows[ dfWindows.G4NN >= dicoParam["g4NN"] ].dropna()
	return dfWindows

def getDFByStrand(df):
	"""Cut the DF with all windows in two, one for each strand.

	:param df: contains all genes windows from one specie.
	:type df: dataFrame

	:returns: dicoDF, contains both reverse and forward dataFrame.
	:rtype: dictionnary of dataFrame
	"""
	dfWindowsReverse = df[ df.geneDesc.str.contains(':-1') ]
	dfWindowsForward = df[ df.geneDesc.str.contains(':[0-9]*:[0-9]*:1',
						regex=True) ]
	dfWindowsReverse = dfWindowsReverse.sort_values(by=['wEnd'])
	dfWindowsForward = dfWindowsForward.sort_values(by=['wStart'])
	dicoDF = {'Reverse' : dfWindowsReverse,
	 		'Forward' : dfWindowsForward}
	return dicoDF

def mergeG4(df, dicoParam, feature, option):
	"""Browses all junction window to find those that are overlapping.

	Here we browse all junctions windows. We will only kept those that overlap
	the 100 nucleotid. Indeed, if the window over thresholds don't overlap this
	position, it only in a gene and not a junction.

	:param df: contains all windows.
	:type df: dataFrame
	:param dicoParam: contains all parameters that were given to g4rna screener.
	:type dicoParam: dictionnary
	:param feature: junction or gene.
	:type feature: string

	:returns: dfpG4, contain all pG4 for that strand.
	:rtype: dataFrame
	"""
	dfTmp = pd.DataFrame()
	dfpG4 = pd.DataFrame()
	dfTmp = dfTmp.append(df[0:1]) # store the first window
	if len(df) == 1:
		dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
				feature, dicoParam["junctionLength"], option))
		dfpG4 = dfpG4.append(dfTmp)
	else:
		for w in range(1,len(df)): # w for window
			# browses all windows over thresholds, exept the first one
			if (df.geneDesc.iloc[w] == df.geneDesc.iloc[w-1] and
			  (feature == 'Gene' and
			        ((df.wStart.iloc[w] >= df.wStart.iloc[w-1] and df.wStart.iloc[w] <= df.wEnd.iloc[w-1]) or
			        (df.wEnd.iloc[w] >= df.wStart.iloc[w-1] and df.wEnd.iloc[w] <= df.wEnd.iloc[w-1])) or
			  feature == 'Junction' and
			        (df.wStart.iloc[w] >= df.wStart.iloc[w-1] and df.wStart.iloc[w] <= df.wEnd.iloc[w-1]))):
				# if window overlap, add window at the current pG4
				dfTmp = dfTmp.append(df[w:w+1])
				if w == len(df)-1 :
					dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
							feature, dicoParam["junctionLength"], option))
					dfpG4 = dfpG4.append(dfTmp)
			else: # new pG4
				dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
						feature, dicoParam["junctionLength"], option))
				dfpG4 = dfpG4.append(dfTmp)
				dfTmp = df.iloc[w:w+1]
				if w == len(df)-1 :
					dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
							feature, dicoParam["junctionLength"], option))
					dfpG4 = dfpG4.append(dfTmp)
				# reinitiate the dfTmp with the new pG4
	return dfpG4

def main(filename, dicoParam, feature, option):
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
		dfpG4 = dfpG4.append(mergeG4(dfWindows, dicoParam, feature, option))
		return dfpG4

if __name__ == '__main__':
	filename = '/home/anais/Documents/Data/Test/pG4.csv'
	dicoParam = {'g4NN' : 0.5, 'cGcC' : 4.5, 'g4H' : 0.9,
				"junctionLength" : 100,
				"windowLength" : 60,
				"step" : 10}
	feature = 'Gene'
	df = main(filename, dicoParam, feature)
	print df.reset_index(drop=True)
