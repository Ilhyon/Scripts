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

def mapG4OnTr(dfTr, pG4):
	"""Finds the location of a pG4r in a transcript.

	Browses all feature of a transcript to find the location of a pG4r. Possible
	locations are : intron, CDS for coding transcript, exon for non coding
	transcript, 3UTR, 5UTR, junction (between two exon), overlap between one
	of these location except junction.

	:param dfTr: contains all feature of a transcript.
	:type dfTr: dataFrame
	:param pG4: contain one pG4r, the one we are checking.
	:type pG4: dataFrame

	:returns: location, location of the pG4 in the transcript.
	:rtype: string
	"""
	location = ' '
	dfTr = dfTr.sort_values(by=['Start'])
	dfTr = dfTr.reset_index(drop=True)
	if len(dfTr) == 1:
		location = dfTr.Feature[0]
	else:
		for index, row in dfTr.iterrows():
			if index != 0:
				if (dfTr.Start[index -1] <= pG4.Start and
					pG4.Start <= dfTr.End[index -1] and
					dfTr.Start[index -1] <= pG4.End and
					pG4.End <= dfTr.End[index -1]):
					#pG4 in the last location
					location = dfTr.Feature[index -1]
					break
				elif (dfTr.End.iloc[index -1] >= pG4.Start and
					pG4.End >= dfTr.Start[index]):
					location = 'overlap_'+ str(dfTr.Feature[index -1]) \
								+'_'+ str(dfTr.Feature[index])
					break
				elif (pG4.Start > dfTr.End[index -1] and
					pG4.End < dfTr.Start[index]):
					location = 'intron'
					break
				elif (pG4.Start <= dfTr.End[index -1] and
					pG4.End > dfTr.End[index -1] and
					pG4.End < dfTr.Start[index]):
					location = 'overlap_'+ str(dfTr.Feature[index -1]) \
								+'_intron'
					break
				elif (pG4.Start > dfTr.End[index -1] and
					pG4.Start < dfTr.Start[index] and
					pG4.End >= dfTr.Start[index]):
					location = 'overlap_intron_' + str(dfTr.Feature.iloc[index])
					break
	return location

def mapG4onJunction(pG4r, dfIntron, dfTr):
	dfpG4Jun = pd.DataFrame()
	for index, row in dfIntron.iterrows():
		dftmp = pG4r
		tr = dfTr[dfTr.Transcript == row.Transcript]
		tr = tr.reset_index(drop=True)
		dftmp['Biotype'] = tr.Biotype[0]
		dfpG4Jun = dfpG4Jun.append(dftmp)
	return dfpG4Jun

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

def main(dfTr, dicoGene, dfpG4, dfIntron):
	dfpG4Tr = pd.DataFrame()
	trRemove = []
	dfpG4Gene =  dfpG4[ dfpG4.Feature == 'Gene' ]
	dfpG4Junction =  dfpG4[ dfpG4.Feature == 'Junction' ]
	for index, row in dfpG4Gene.iterrows():
		geneId = row.id
		dfGtfGene =  dfTr[ dfTr.Gene == geneId ]
		transcripts = list(set(dfGtfGene.Transcript))
		for tr in transcripts:
			dfGtfTr =  dfGtfGene[ dfGtfGene.Transcript == tr ].dropna()
			if row.Start >= min(dfGtfTr.Start) and row.End <= max(dfGtfTr.End):
				# the pG4 is in this transcript
				if dfGtfTr.Type.iloc[0] == 'Coding':
					dfGtfTr = dfGtfTr[dfGtfTr.Feature != 'exon']
				else:
					if len(dfGtfTr[ dfGtfTr.Feature.str.contains('utr') ]) > 0:
						# tr with bad annotation
						trRemove.append(tr)
				location = mapG4OnTr(dfGtfTr, row)
				pG4rtmp = row
				pG4rtmp['Location'] = location
				pG4rtmp['Biotype'] = list(set(dfGtfTr.Biotype))[0]
				pG4rtmp.id = tr
				dfpG4Tr = dfpG4Tr.append(pG4rtmp)
	pG4rtmp = pd.DataFrame()
	for index, row in dfpG4Junction.iterrows():
		id = row.id
		dftmpIntron = dfIntron[ dfIntron.Id == id ]
		dftmpIntron = dftmpIntron.reset_index(drop=True)
		pG4rtmp = pG4rtmp.append(mapG4onJunction(dfpG4Junction, dftmpIntron, dfTr))
	if len(pG4rtmp) > 0:
		pG4rtmp = pG4rtmp.drop_duplicates(subset=None, keep='first', inplace=False)
		print(pG4rtmp)
		pG4rtmp['Location'] = 'junction'
		dfpG4Tr = dfpG4Tr.append(pG4rtmp)
		dfpG4Tr = dfpG4Tr.reset_index(drop=True)
	dfpG4Tr = removeG4OnBadTr(dfpG4Tr, trRemove)
	dfpG4Tr = dfpG4Tr.reset_index(drop=True)
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
