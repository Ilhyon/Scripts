#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
import Parser_gtf
import numpy as np
import pandas as pd

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	# GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = '/home/anais/Documents/Data/Genomes')
	parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
	return parser

def overlaps(interval1, interval2):
    """Computes the distance or overlap between two interval.

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
	"""Gets the location of the pG4 in a transcript.

	We compute the overlap between an exon and a pG4. If there is an overlap,
	we search if the pG4 is overlaping an an intron or if it is entirely in
	the exon. If not, it's in the itron.

	:param coordExon: coords of one exon.
	:type coordExon: list
	:param coordpG4: coords of a pG4.
	:type coordpG4: list

	:returns: location, location of the pG4
	:rtype: string
	"""
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
		location = None
	return location

def mapOnTr(dfpG4gene, dfGTF, dicoGene):
	"""Maps all pG4 that are in gene on transcript and location.

	All gene's transcript that contain pG4 are browse to find the pG4 location.
	There is 5 type of location : overlap Intron-Exon, Intron, Exon and
	multi overlap if is on more thant 2 feature (example : F01D5.10.1 with the
	pG4 [14019141, 14019610]).

	:param dfpG4gene: conatin all gene's pG4.
	:type dfpG4gene: dataFrame
	:param dfGTF: GTF file informations.
	:type dfGTF: dataFrame


	:returns: dfpG4Tr, all pG4 at the transcript level, with their location.
	:rtype: dataFrame
	"""
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
				location = [mapPremRNA(coordExon, coordpG4) for coordExon in dfTrexon.Coords]
				if len(list(set(location))) == 1 and location[0] == None:
					location = 'Intron'
				elif len(list(set(location))) > 1:
					location = 'Multi_overlap'
				else:
					location = list(set([ loc for loc in location if loc != None]))
				dfTmp = pd.DataFrame.from_dict({'Transcript' : tr,
						'Location' : location, 'Sequence' : row.seqG4,
						'cGcC' : row.cGcC, 'G4H' : row.G4H, 'G4NN' : row.G4NN,
						'Biotype' : biotype, 'Type' : type})
				dfpG4Tr = dfpG4Tr.append(dfTmp)
	return dfpG4Tr

def mapOnJunction(dfpG4Junction, dfGTF, dfIntron):
	"""Maps all pG4 that are in junctions on transcript.

	Initialy pG4 on junction are uniq, because they are only identified by their
	chromosomal coords. So here don't find their location but only their
	transcripts.

	:param dfpG4Junction: conatin all junctions's pG4.
	:type dfpG4Junction: dataFrame
	:param dfGTF: GTF file informations.
	:type dfGTF: dataFrame
	:param dfIntron: contain all intron and their informations.
	:type dfIntron: dataFrame

	:returns: dfpG4Tr, all pG4 at the transcript level, with their location.
	:rtype: dataFrame
	"""
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
	"""Finds location and transcript of G4 predicted in genes.

	We starts by pG4 in genes and then in junction. They are only mapped at the
	pre-mRNA level (Exon-Intron) because a lot species do not pocess the UTR
	annotation.

	:param dfGTF: GTF file informations.
	:type dfGTF: dataFrame
	:param dicoGene: contain all gene and their transcripts {idGene : idTr}
	:type dicoGene: dictionary
	:param dfpG4: conatin all predicted G4 at the gene level.
	:type dfpG4: dataFrame
	:param dfIntron: contain all intron and their informations.
	:type dfIntron: dataFrame
	"""
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
