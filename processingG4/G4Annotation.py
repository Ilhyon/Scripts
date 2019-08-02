#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
import Parser_gtf
import numpy as np
import pandas as pd
import timeit as ti
from pprint import pprint

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	# GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = '/home/anais/Documents/Data/Genomes')
	parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
	return parser

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

def mapPremRNA(coordExon, coordpG4):
	o = overlaps(coordExon, coordpG4)
	if o > 0:
		# the pG4 overlap the exon
		if o >= (coordpG4[1] - coordpG4[0]):
			# the entire pG4 is overlaping the exon
			location = 'Exon'
		elif coordpG4[0] < coordExon[0]:
			location = 'Overlap_Intron_Exon'
		else:
			location = 'Overlap_Exon_Intron'
	else:
		location = 'Intron'
	return location

def removeG4OnBadTr(dfpG4Tr, trRemove):
	""" Remove from a dataFrame row given by a list.

	Some transcript get a bad annotation (UTR in non coding), so we need to
	remove them from our data.

	:param dfpG4Tr: all pG4r annotated.
	:type dfpG4Tr: dataFrame
	:param trRemove: all transcript we need to remove.
	:type trRemove: list

	:returns: dfpG4Tr, filtered dataFrame with all pG4r.
	:rtype: dataFrame
	"""
	for index, row in dfpG4Tr.iterrows():
		tr = row.id
		if tr in trRemove:
			dfpG4Tr = dfpG4Tr[dfpG4Tr.id != tr]
	return dfpG4Tr

def mapOnTr(dfpG4gene, dfGTF, dicoGene):
	dfpG4Tr = pd.DataFrame()
	# retrieve gene with G4 and then there transcripts
	geneWithpG4 = list(set(dfpG4gene.id))
	transcripts = [tr for g in geneWithpG4 for tr in dicoGene[g]]
	for tr in transcripts:
		# retrieve tr exons
		dfTr =  dfGTF[ dfGTF.Transcript == tr ].dropna()
		dfTrexon = dfTr[ dfTr.Feature == 'exon'].dropna()
		gene = list(set(dfTr.Gene))[0]
		biotype =  list(set(dfTrexon.Biotype))
		type = list(set(dfTrexon.Type))
		# retrieve all pG4 in the current gene
		pG4inGene = dfpG4gene[ dfpG4gene.id == gene]
		for index, row in pG4inGene.iterrows():
			if row.Start >= min(dfTr.Start) and row.End <= max(dfTr.End):
				# the current pG4 is in the transcripts
				coordpG4 = [row.Start, row.End]
				location = [mapPremRNA(coordExon, coordpG4) for coordExon in dfTrexon.Coords][0]
				dfTmp = pd.DataFrame.from_dict({'Transcript' : tr,
						'Location' : location, 'Sequence' : row.seqG4,
						'cGcC' : row.cGcC, 'G4H' : row.G4H, 'G4NN' : row.G4NN,
						'Biotype' : biotype, 'Type' : type})
				dfpG4Tr = dfpG4Tr.append(dfTmp)
	return dfpG4Tr

def mapOnJunction(dfpG4Junction, dfGTF, dfIntron):
	dfpG4Tr = pd.DataFrame()
	# retrieve introns that have pG4
	introns = list(set(dfpG4Junction.id))
	for i in introns:
		# get all transcript that get this intron (thus the junction E-E)
		transcripts = list(set(dfIntron[dfIntron.Id == i].Transcript))
		for tr in transcripts:
			biotype = list(set(dfGTF[ dfGTF.Transcript == tr].Biotype))[0]
			type = list(set(dfGTF[ dfGTF.Transcript == tr].Type))[0]
			pG4 = dfpG4Junction[ dfpG4Junction.id == i]
			dfTmp = pd.DataFrame.from_dict({'Transcript' : tr,
					'Location' : 'Junction', 'Sequence' : pG4.seqG4,
					'cGcC' : pG4.cGcC, 'G4H' : pG4.G4H, 'G4NN' : pG4.G4NN,
					'Biotype' : biotype, 'Type' : type})
			dfpG4Tr = dfpG4Tr.append(dfTmp)
	return dfpG4Tr


def main(dfGTF, dicoGene, dfpG4, dfIntron):
	dfpG4Tr = pd.DataFrame()
	# split G4 in gene from G4 in junction
	dfpG4gene =  dfpG4[ dfpG4.Feature == 'Gene' ]
	dfpG4Junction =  dfpG4[ dfpG4.Feature == 'Junction' ]
	dfpG4Tr = dfpG4Tr.append(mapOnTr(dfpG4gene, dfGTF, dicoGene))
	dfpG4Tr = dfpG4Tr.append(mapOnJunction(dfpG4Junction, dfGTF, dfIntron))
	return dfpG4Tr

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	sp = arg.specie
	filename = path+'/'+sp+'/'+sp+'.gtf'
	dfTr = Parser_gtf.importGTFdf(filename)
	dicoGene = Parser_gtf.importGTFGene(filename)
	main(dfTr, dicoGene, dfpG4)
