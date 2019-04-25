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

def getChromosomalPositionForwardStrand(position,
									junLength,
									startIntron,
									endIntron):
	"""
		Give chromosomal positions of G4 around a junction,
		in a forward gene.
	"""
	if position <= junLength:	# upstream sequence
		upstreamLength = junLength - position
		position = startIntron - 1 - upstreamLength
	else:	# downstream sequence
		downstreamLength = position - junLength
		position = endIntron + 1 + downstreamLength - 1
	return position

def getChromosomalPositionReverseStrand(position,
									junLength,
									startIntron,
									endIntron):
	"""
		Give chromosomal positions of G4 around a junction,
		in a reverse gene.
	"""
	if position <= junLength + 1:	# because strart G4 classifier from 0
		position = startIntron + 1 + junLength + 1 - position
	else:	# if from sequence aval
		position = endIntron + junLength - position
	return position

def getChromosomalCoord(pG4Start, pG4End, geneDesc, junctionLength):
	"""
		This function is principaly used for junction.
	"""
	intronStart = int(geneDesc.split("|")[1].split("-")[0])
	intronEnd = int(geneDesc.split("|")[1].split("-")[1])
	if intronStart < intronEnd: # gene on forward strand
		pG4Start = getChromosomalPositionForwardStrand(pG4Start,
					junctionLength,
					intronStart,
					intronEnd)
		pG4End = getChromosomalPositionForwardStrand(pG4Start,
					junctionLength,
					intronStart,
					intronEnd)

	else: # gene on reverse strand
		pG4Start = getChromosomalPositionReverseStrand(pG4Start,
					junctionLength,
					intronStart,
					intronEnd)
		pG4End = getChromosomalPositionReverseStrand(pG4Start,
					junctionLength,
					intronStart,
					intronEnd)
		tmp = pG4Start
		pG4Start = pG4End
		pG4End = tmp
	return pG4Start, pG4End

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

def mergeWindows(dfTmp, feature, junctionLength):
	"""Merge overlaping windows.

	:param dfTmp: contain overlaping windows.
	:type dfTmp: dataFrame
	:param feature: junction or gene.
	:type feature: string
	:param junctionLength: length of unction, by default it's 100 nt.
	:type junctionLength: integer

	:returns: pG4, contains the pG4 which is the merge of overlaping windows.
	:rtype: dictionary
	"""
	lastRow = len(dfTmp.index) - 1
	geneDesc = dfTmp.geneDesc.iloc[0]
	geneId = geneDesc.split(' ')[0]
	geneDescSplit = geneDesc.split(' ')[2].split(':')
	strand = geneDescSplit[5]
	geneStart = int(geneDescSplit[3])
	geneEnd = int(geneDescSplit[4])
	meancGcC = dfTmp.cGcC.mean()
	meanG4H = dfTmp.G4H.mean()
	meanG4NN = dfTmp.G4NN.mean()
	pG4Start = int(dfTmp.wStart.iloc[0])
	pG4End = int(dfTmp.wEnd.iloc[lastRow])
	pG4rSeq = mergeOverlappingSequences(dfTmp)
	if feature == "Junction":
		if pG4Start < junctionLength and pG4End > junctionLength :
			pG4Start, pG4End = getChromosomalCoord(pG4Start, pG4End,
								geneDesc, junctionLength)
			pG4 = {"geneDesc" : [geneDesc], "Gene" : geneId,
					"cGcC" : [meancGcC], "G4H" : [meanG4H], "G4NN" : [meanG4NN],
					"pG4Start" : [pG4Start], "pG4End" : [pG4End],
					"seqG4" : [pG4Seq],
					"Feature" : [feature]}
		else:
			pG4 = None
	else:
		pG4 = {"geneDesc" : [geneDesc], "Gene" : geneId,
				"cGcC" : [meancGcC], "G4H" : [meanG4H], "G4NN" : [meanG4NN],
				"pG4Start" : [pG4Start], "pG4End" : [pG4End],
				"seqG4" : [pG4rSeq],
				"Feature" : [feature]}
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
	dfWindowsForward = df[ df.geneDesc.str.contains(':[0-9]*:[0-9]*:1', regex=True) ]
	dfWindowsReverse = dfWindowsReverse.sort_values(by=['wStart'])
	dfWindowsForward = dfWindowsForward.sort_values(by=['wEnd'])
	dicoDF = {'Reverse' : dfWindowsReverse,
	 		'Forward' : dfWindowsForward}
	return dicoDF

def merge(df, dicoParam, feature):
	"""Browses all window to find those that are overlapping.

	:param df: contains all windows from one strand.
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
	if len(df) == 1 :
		dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
				feature, dicoParam["junctionLength"]))
		dfpG4 = dfpG4.append(dfTmp)
	else:
		for w in range(1,len(df)): # w for window
			# browses all windows over thresholds, exept the first one
			if ((df.wStart.iloc[w] >= df.wStart.iloc[w-1] and
				df.wStart.iloc[w] <= df.wEnd.iloc[w-1]) or
				(df.wEnd.iloc[w] >= df.wStart.iloc[w-1] and
				df.wEnd.iloc[w] <= df.wEnd.iloc[w-1])):
				# if window overlap, add window at the current pG4
				dfTmp = dfTmp.append(df[w:w+1])
				if w == len(df)-1 :
					dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
							feature, dicoParam["junctionLength"]))
					dfpG4 = dfpG4.append(dfTmp)
			else: # new pG4
				dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
						feature, dicoParam["junctionLength"]))
				dfpG4 = dfpG4.append(dfTmp)
				dfTmp = df.iloc[w:w+1]
				if w == len(df)-1 :
					dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
							feature, dicoParam["junctionLength"]))
					dfpG4 = dfpG4.append(dfTmp)
				# reinitiate the dfTmp with the new pG4
	return dfpG4

def main(filename, dicoParam, feature):
	dfpG4 = pd.DataFrame()
	try:
		dfWindows = pd.read_csv(filename, sep='\t', index_col=0)
	except:
		print 'Buuuuuuug ' +filename
	else:
		# dataFrame with all windows from G4RNA Screener
		dfWindows.columns = ['geneDesc','cGcC',
							'G4H','seqG4','wStart',
							'wEnd', 'G4NN']
		dfWindows = filterOnScores(dicoParam, dfWindows)
		dfWStrand = getDFByStrand(dfWindows)
		for strand in dfWStrand:
			dfpG4 = dfpG4.append(merge(dfWStrand[strand], dicoParam, feature))
		return dfpG4

if __name__ == '__main__':
	main(filename, dicoParam, feature)
