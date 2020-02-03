#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	August 2019

Description:
	This file is more like a library called when statistics are computed (
		getMainDensities and getDataFig). Its aims is to compute the number of
		location with at least one pG4 in it. A special script is made for
		because we need to compute each time, otherwise we could get the same
		location in different condition and count it 2 times.
"""

import os
import argparse
import getMainDensities
from pprint import pprint
import recurrentFunction as rF

def overlaps(interval1, interval2):
	"""Compute the distance or overlap between two interval.

	Compute : min(ends) - max(starts). If > 0, return the number of bp of
	overlap, if 0, they are book-ended and if < 0 return the distance in
	bp between them.

	:param interval1: contain the start and the end of a OQs.
	:type interval1: list of int
	:param interval2: contain the start and the end of a gff feature.
	:type interval2: list of int

	:returns: min(ends) - max(starts)
	:rtype: int
	"""
	return min(interval1[1], interval2[1]) - max(interval1[0], interval2[0])

def overlapLocation(coordspG4, coordsLoc):
	"""Gets the overlap between two set of coordinates.

	If there is an overlap (>0) between an pG4 and the location, then the
	coords of the location are returned.

	:param intervalSeq: contain the start and the end of a pG4.
	:type intervalSeq: list of int
	:param intervalFeat: contain the start and the end of a gff feature.
	:type intervalFeat: list of int

	:returns: ':'.join(coordsLoc), location coord if overlap with the pG4.
	:rtype: list
	"""
	#if element are on the same chr and strand
	coordspG4 = coordspG4.split(':')
	coordsLoc = coordsLoc.split(':')
	if coordspG4[0] == coordsLoc[0] and coordspG4[2] == coordsLoc[2]:
		coordG4 = coordspG4[1].split('~')
		coordG4 = [int(coordG4[0]), int(coordG4[1]) ]
		coordLoc = coordsLoc[1].split('~')
		coordLoc = [int(coordLoc[0]), int(coordLoc[1]) ]
		overlap = overlaps(coordG4, coordLoc)
		if overlap > 0:
			return ':'.join(coordsLoc)

def importLocaMulti(filename):
	"""Import into a dictionary the number of uniq location by transcript.

	pG4 are affiliated to a chromosomal location and a transcript. So here we
	count the number tot of location according to transcripts. To do that all
	locations are browsed, then transcript

	:param filename: Location filename, this file contain all locations
		coordinates and all transcript-biotype that 'possess' it.
	:type filename: string

	:returns: dicoLoca, {Location-Biotype : [idLocation-idTr]}
	:rtype: dictionary
	"""
	dicoLoca = {}
	with open(filename) as f:
		lines = f.read().splitlines()
		for l in lines:
			l = l.rstrip()
			w = l.split(':')
			w[0] = w[0].split('>')[1]
			if w[5]:
				location = w[1]
				listTrBt = w[5].split('|')
				chr = w[2]
				coords = w[3]
				strand = w[4]
				dicoTrBt = { TrBt.split('-')[0] : TrBt.split('-')[1] for TrBt in listTrBt}
				idLoca = chr+':'+coords+':'+strand
				for tr in dicoTrBt:
					locBt = location+'-'+dicoTrBt[tr]
					idLoca = idLoca+tr
					if locBt not in dicoLoca:
						dicoLoca[locBt] = []
					dicoLoca[locBt].append(idLoca)
	return dicoLoca

def importLoca(filename):
	"""Import into a dictionary the number of uniq location.

	Same as importLocaMulti but without the transcript constrain.

	:param filename: Location filename, this file contain all locations
		coordinates and all transcript-biotype that 'possess' it.
	:type filename: string

	:returns: dicoLoca, {Location-Biotype : [idLocation]}
	:rtype: dictionary
	"""
	dicoLoca = {}
	with open(filename) as f:
		lines = f.read().splitlines()
		for l in lines:
			l = l.rstrip()
			w = l.split(':')
			w[0] = w[0].split('>')[1]
			if w[5]:
				location = w[1]
				listTrBt = w[5].split('|')
				chr = w[2]
				coords = w[3]
				strand = w[4]
				listBt = [ TrBt.split('-')[1] for TrBt in listTrBt ]
				idLoca = chr+':'+coords+':'+strand
				for bt in listBt:
					locBt = location+'-'+bt
					if locBt not in dicoLoca:
						dicoLoca[locBt] = []
					dicoLoca[locBt].append(idLoca)
	return dicoLoca

def countLocId(dicoLocId):
	"""Gets the number of location from a dictionary {loc-Bt : [locID]}

	This function just browse the input dictionary and for each key (loc-Bt),
	make the list unique and its length is the number of locations.

	:param dicoLocId: {Location-Biotype : [locId]}
	:type dicoLocId: dictionary
	"""
	dicoLoca = {}
	for locBt in dicoLocId:
		dicoLoca[locBt] = len(set(dicoLocId[locBt]))
	return dicoLoca

def countLocIdByClass(dicoLocId):
	"""Get the number of location by class from biotype.

	The dictionary with all locId is parser to another one where keys are
	changed from the biotype level to the class level (coding, long ncRNA, short
	ncRNA, and pseudogene).

	:param dicoLocId: biotype level of location list
	:type dicoLocId: dictionary

	:returns: dicoLoca, class level.
	:rtype: dictionary
	"""
	dicoLoca = {}
	dicoClass = rF.createDicoFamilyFiltered()
	for locBt in dicoLocId:
		loc = locBt.split('-')[0]
		Bt = locBt.split('-')[1]
		for classe in dicoClass:
			if Bt in dicoClass[classe]:
				curClass = classe
				locClass = loc+'-'+curClass
				if locClass not in dicoLoca:
					dicoLoca[locClass] = []
				dicoLoca[locClass].extend(dicoLocId[locBt])
	for locClass in dicoLoca:
		dicoLoca[locClass] = len(set(dicoLoca[locClass]))
	return dicoLoca

def readShufpG4(pG4File):
	"""Reads a file containing all pG4 Shuf and parse it.

	This function aims to retrieve for each biotype location, the location where
	a pG4 is.

	:param filename: name of the file containing all shuf pG4.
	:type filename: string

	:returns: dicoLoca, {Location-Biotype : [coordLocation]}
	:rtype: dictionary
	"""
	dicoLoca = {}
	with open(pG4File) as f:
		lines = f.read().splitlines()
		for l in lines:
			l = l.rstrip()
			w = l.split('\t')
			if w[0] != 'Transcript' and w:
				location = w[1]
				tr = w[0]
				bt = w[10]
				chr = w[4]
				start = w[3]
				end = w[5]
				strand = w[6]
				idLoca = chr + ':' + start + '~' + end + ':' + strand
				locBt = location+'-'+bt
				if locBt not in dicoLoca:
					dicoLoca[locBt] = []
				dicoLoca[locBt].append(idLoca)
	return dicoLoca

def readWtpG4(filename, dicoLoc):
	"""Reads the pG4 Wt file and map them on location to find their coords.

	In the input file, locations name are given but not their coords. So we need
	to map them to make them unique. Each line are red, and for each pG4,
	we find the overlap between the pG4 and all location of this pG4 in this
	transcript. For example, there is a pG4 in a transcript 1. This Tr 1 have
	4CDS, 1 5UTR, 1 3UTR and 3 introns. So if this pG4 is annotated in a CDS,
	then all CDS will be browsed to find in wich one is the pG4. Then this CDS
	and its coordinates are saved in a list.
	In this file, exon coding are not registered. So when a pG4 in a UTR, CDS
	or codon is found, the exon is also 'created' to don't loose this
	information.

	:param filename: name of the file with all Wt pG4 annotated.
	:type filename: string
	:param dicoLoc: dicoLoc, {idTr : {LocationName : [coords]} }
	:type dicoLoc: dictionary

	:returns: dicoLoca, {Location-Biotype : [coordLocation]}
	:rtype: dictionary
	"""
	dicoLoca = {}
	with open(filename) as f:
		lines = f.read().splitlines()
		for l in lines:
			l = l.rstrip()
			w = l.split('\t')
			if w[0] != 'Transcript':
				start = w[3]
				end = w[5]
				chr = w[4]
				strand = w[6]
				coordspG4 =  chr +':'+ start +'~'+ end +':'+ strand
				loc = w[1]
				trId = w[0]
				locBt = loc+'-'+w[10]
				if loc in dicoLoc[trId]:
					overlaps = [ overlapLocation(coordspG4, coordsLoc) for coordsLoc in dicoLoc[trId][loc] ]
					listCommon = [x for x in overlaps if x is not None]
					if locBt not in dicoLoca:
						dicoLoca[locBt] = []
					dicoLoca[locBt].extend(listCommon)
				if loc in ['5UTR', 'CDS', '3UTR', 'StartCodon', 'StopCodon']:
					locBt = 'exon-'+w[6]
					loc = 'exon'
					overlaps = [ overlapLocation(coordspG4, coordsLoc) for coordsLoc in dicoLoc[trId][loc] ]
					listCommon = [x for x in overlaps if x is not None]
					if locBt not in dicoLoca:
						dicoLoca[locBt] = []
					dicoLoca[locBt].extend(listCommon)
	return dicoLoca

def importIndexByTr(filename):
	"""Creates a dictionary with all Tr and all their location.

	Browses the input location file and create a dictionary of dictionary :
	{idTr : {LocationName : [coords]} }. So we have for all transcript, all
	their location, and for all location we have their coordinates.

	:param filename: Location filename, this file contain all locations
		coordinates and all transcript-biotype that 'possess' it.
	:type filename: string

	:returns: dicoLoca, {idTr : {LocationName : [coords]} }.
	:rtype: dictionary
	"""
	dicoLoca = {}
	with open(filename) as f:
		lines = f.read().splitlines()
		for l in lines:
			l = l.rstrip()
			w = l.split(':')
			w[0] = w[0].split('>')[1]
			if w[5]:
				location = w[1]
				listTrBt = w[5].split('|')
				chr = w[2]
				coords = w[3]
				strand = w[4]
				listTr = [ TrBt.split('-')[0] for TrBt in listTrBt ]
				for tr in listTr:
					if tr not in dicoLoca:
						dicoLoca[tr] = {}
					if location not in dicoLoca[tr]:
						dicoLoca[tr][location] = []
					dicoLoca[tr][location].append(chr+':'+coords+':'+strand)
	return dicoLoca

def getAllLocNum(path, type):
	"""Creates a dictionary with the number of location.

	Same as getAllLocNumMulti but this time location are not linked to
	transcripts. First data are imported into a dictionary
	{Location-Biotype : [Location-Tr]}. Then, depending on if we want the
	result at the class or biotype level, the list of location is made
	unique and its length give the number of location.
	This result is also store in a dictionary {LevelLocation : numberOfLocation}.

	:param path: main folder path with all needed file in it.
	:type path: string
	:param type: level (class or biotype/subclass)
	:type type: string

	:returns: dicoNbTot, {LevelLocation : numberOfLocation}
	:rtype: dictionary
	"""
	fileLoca = path+'/Locations.txt'
	dicoLocID = importLoca(fileLoca)
	if type == 'Class':
		dicoNbTot = countLocIdByClass(dicoLocID)
	else:
		dicoNbTot = countLocId(dicoLocID)
	return dicoNbTot

def getAllLocNumMulti(path, type):
	"""Creates a dictionary with the number of location.

	pG4 are affiliated to a chromosomal location and a transcript. So here we
	count the number tot of location according to transcripts. First data are
	imported into a dictionary {Location-Biotype : [Location-Tr]}. Then,
	depending on if we want the result at the class or biotype level, the
	list of location is made unique and its length give the number of location.
	This result is also store in a dictionary {LevelLocation : numberOfLocation}.

	:param path: main folder path with all needed file in it.
	:type path: string
	:param type: level (class or biotype/subclass)
	:type type: string

	:returns: dicoNbTot, {LevelLocation : numberOfLocation}
	:rtype: dictionary
	"""
	fileLoca = path+'/Locations.txt'
	dicoLocID = importLocaMulti(fileLoca)
	if type == 'Class':
		dicoNbTot = countLocIdByClass(dicoLocID)
	else:
		dicoNbTot = countLocId(dicoLocID)
	return dicoNbTot

def getpG4WtLocNum(path, type):
	"""Creates a dictionary with the number of location with at least one pG4.

	This function aim to find the number of location with at least one pG4 in
	it.

	:param path: main folder path with all needed file in it.
	:type path: string
	:param type: level (class or biotype/subclass)
	:type type: string

	:returns: dicoNbTot, {LevelLocation : numberOfLocation}
	:rtype: dictionary
	"""
	fileLoca = path+'/Locations.txt'
	filepG4Wt = path+'/pG4_shuffled.csv'
	dicoLocTr = importIndexByTr(fileLoca)
	dicopG4Wt = readWtpG4(filepG4Wt, dicoLocTr)
	if type == 'Class':
		dicoNbWt = countLocIdByClass(dicopG4Wt)
	else:
		dicoNbWt = countLocId(dicopG4Wt)
	return dicoNbWt

def getpG4ShufLocNum(path, type):
	"""Gets the number of location with at least one pG4 in the Shuf dataset.

	This function aim to find the number of location with at least one pG4 in
	it.

	:param path: main folder path with all needed file in it.
	:type path: string
	:param type: level (class or biotype/subclass)
	:type type: string

	:returns: dicoNbTot, {LevelLocation : numberOfLocation}
	:rtype: dictionary
	"""
	filepG4Shuf = path +'/pG4_shuffled.csv'
	dicopG4Shuf = readShufpG4(filepG4Shuf)
	if type == 'Class':
		dicoNbShuf = countLocIdByClass(dicopG4Shuf)
	else:
		dicoNbShuf = countLocId(dicopG4Shuf)
	return dicoNbShuf

def main(path):
	"""Main function, only used for tests.
	"""
	fileLoca = mainPath+'/Locations.txt'
	filepG4Wt = mainPath+'/HS_All_G4InTranscript.txt'
	filepG4Shuf = mainPath+'/pG4_shuffle.csv'
	pprint(getAllLocNumMulti(path, 'bt'))
	# dicoLocID = importLoca(fileLoca)
	# pprint(countLocIdByClass(dicoLocID))
	# pprint(countLocId(dicoLocID))

	# dicoLocTr = importIndexByTr(fileLoca)
	# pprint(countLocIdByClass(dicoLocTr))
	# dicopG4Wt = readWtpG4(filepG4Wt, dicoLocTr)
	# pprint(countLocId(dicopG4Wt))
	# dicopG4Shuf = readShufpG4(filepG4Shuf)
	# pprint(countLocId(dicopG4Shuf))

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'analyseGC')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	mainPath = arg.path
	main(mainPath)
