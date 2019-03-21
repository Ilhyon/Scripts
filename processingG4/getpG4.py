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

def mergeOverlappingSequences(dfTmp):
	seq = str(dfTmp.seqG4.iloc[0])
	for w in range(1,len(dfTmp)):
		step = int(dfTmp.wStart.iloc[w] - dfTmp.wStart.iloc[w-1])
		# convert to int elsewise it's a float
		wSeq = dfTmp.seqG4.iloc[w]
		seq += wSeq[-step:]
	return seq

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

def changeCoordByStrand(start, end):
	if start > end: # gene on reverse strand
		tmp = start
		start = end
		end = tmp
	return start, end

def mergeWindows(dfTmp, feature, junctionLength):
	"""
		Merge all widnows from a dataFrame to a unique pG4
	"""
	lastRow = len(dfTmp.index) -1
	geneDesc = dfTmp.geneDesc.iloc[0]
	geneId = geneDesc.split("|")[0]
	meancGcC = dfTmp.cGcC.mean()
	meanG4H = dfTmp.G4H.mean()
	meanG4NN = dfTmp.G4NN.mean()
	pG4Start = dfTmp.wStart.iloc[0]
	pG4End = dfTmp.wEnd.iloc[lastRow]
	pG4Start, pG4End = changeCoordByStrand(pG4Start, pG4End)
	if feature == "Junction":
		if pG4Start < junctionLength and pG4End > junctionLength :
			pG4Start, pG4End = getChromosomalCoord(pG4Start, pG4End,
								geneDesc, junctionLength)
			pG4Seq = mergeOverlappingSequences(dfTmp)
			pG4 = {"geneDesc" : [geneDesc], "Gene" : geneId,
					"cGcC" : [meancGcC], "G4H" : [meanG4H], "G4NN" : [meanG4NN],
					"pG4Start" : [pG4Start], "pG4End" : [pG4End],
					"seqG4" : [pG4Seq],
					"Feature" : [feature]}
		else:
			pG4 = None
	else:
		pG4Seq = mergeOverlappingSequences(dfTmp)
		pG4 = {"geneDesc" : [geneDesc], "Gene" : geneId,
				"cGcC" : [meancGcC], "G4H" : [meanG4H], "G4NN" : [meanG4NN],
				"pG4Start" : [pG4Start], "pG4End" : [pG4End],
				"seqG4" : [pG4Seq],
				"Feature" : [feature]}
	return pG4

def filterOnScores(dicoParam, dfWindows):
	"""
		Filters the data frame with all windows, we only keep windiows
		over all thresholds.
	"""
	dfWindows = dfWindows[ dfWindows.cGcC >= dicoParam["cGcC"] ].dropna()
	dfWindows = dfWindows[ dfWindows.G4H >= dicoParam["g4H"] ].dropna()
	dfWindows = dfWindows[ dfWindows.G4NN >= dicoParam["g4NN"] ].dropna()
	return dfWindows

def mainDetectpG4(filename, dicoParam, feature):
	dfpG4 = pd.DataFrame()
	dfTmp = pd.DataFrame()
	dfWindows = pd.read_csv(filename, sep='\t', index_col=0)
	# dataFrame with all windows from G4RNA Screener
	dfWindows.columns = ['geneDesc','cGcC',
						'G4H','seqG4','wStart',
						'wEnd', 'G4NN']
	dfWindows = filterOnScores(dicoParam, dfWindows)
	dfTmp = dfTmp.append(dfWindows[0:1]) # store the first window
	for w in range(1,len(dfWindows)): # w for window
		# ~ print dfWindows[w:w+1] #-> ligne w
		# browses all windows over thresholds, exept the first
		if (dfWindows.wStart.iloc[w] >= dfWindows.wStart.iloc[w-1] and
		dfWindows.wStart.iloc[w] <= dfWindows.wEnd.iloc[w-1] and
		dfWindows.geneDesc.iloc[w] == dfWindows.geneDesc.iloc[w-1]):
			# add window for current pG4
			dfTmp = dfTmp.append(dfWindows[w:w+1])
		else: # new pG4
			dfTmp = pd.DataFrame.from_dict(mergeWindows(dfTmp,
					feature, dicoParam["junctionLength"]))
			dfpG4 = dfpG4.append(dfTmp)
			dfTmp = dfWindows.iloc[w:w+1] 
			# reinitiate the dfTmp with the new pG4
	return dfpG4

if __name__ == '__main__':
	mainDetectpG4(filename, dicoParam, feature)
