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
	if len(dfTr) == 1:
		location = dfTr.Feature[0]
	else:
		for feature in range(1,len(dfTr)):
			if (dfTr.Start.iloc[feature -1] <= pG4.Start[0] and
				pG4.Start[0] <= dfTr.End.iloc[feature -1] and
				dfTr.Start.iloc[feature -1] <= pG4.End[0] and
				pG4.End[0] <= dfTr.End.iloc[feature -1]):
				#pG4 in the last location
				location = dfTr.Feature.iloc[feature -1]
				break
			elif (dfTr.End.iloc[feature -1] >= pG4.Start[0] and
				pG4.End[0] >= dfTr.Start.iloc[feature]):
				location = 'overlap_'+ str(dfTr.Feature.iloc[feature -1]) \
							+'_'+ str(dfTr.Feature.iloc[feature])
				break
			elif (pG4.Start[0] > dfTr.End.iloc[feature -1] and
				pG4.End[0] < dfTr.Start.iloc[feature]):
				location = 'intron'
				break
			elif (pG4.Start[0] <= dfTr.End.iloc[feature -1] and
				pG4.End[0] > dfTr.End.iloc[feature -1] and
				pG4.End[0] < dfTr.Start.iloc[feature]):
				location = 'overlap_'+ str(dfTr.Feature.iloc[feature -1]) \
							+'_intron'
				break
			elif (pG4.Start[0] > dfTr.End.iloc[feature -1] and
				pG4.Start[0] < dfTr.Start.iloc[feature] and
				pG4.End[0] >= dfTr.Start.iloc[feature]):
				location = 'overlap_intron_' + str(dfTr.Feature.iloc[feature])
				break
	return location

def mapG4onJunction(pG4r, dfIntron, dfTr):
	dfpG4Jun = pd.DataFrame()
	for intron in range(0,len(dfIntron)):
		dftmp = pG4r
		dftmp['Location'] = 'junction'
		tr = dfTr[ dfTr.Transcript == dfIntron.Transcript[intron] ]
		tr = tr.reset_index(drop=True)
		dftmp['Biotype'] = tr.Biotype[0]
		dfpG4Jun = dfpG4Jun.append(dftmp)
	return dfpG4Jun

def removeG4OnBadTr(dfpG4Annotation, trRemove):
	for pG4 in range(0,len(dfpG4Annotation)):
		print pG4
		tr = dfpG4Annotation.id[pG4]
		if tr in trRemove:
			dfpG4Annotation = dfpG4Annotation[dfpG4Annotation.id != tr]
	return dfpG4Annotation

def main(dfTr, dicoGene, dfpG4, dfIntron):
	dfpG4Annotation = pd.DataFrame()
	trRemove = []
	for pG4 in range(0,len(dfpG4)):
		geneId = dfpG4.id.iloc[pG4]
		if dfpG4.Feature[pG4] == 'Gene':
			if geneId in dicoGene:
				if 'Transcript' in dicoGene[geneId]:
					for tr in dicoGene[geneId]['Transcript']:
						pG4r = dfpG4.iloc[pG4:pG4+1]
						pG4r = pG4r.reset_index(drop=True)
						dftmp =  dfTr[ dfTr.Transcript == tr ].dropna()
						if dftmp.Type.iloc[0] == 'Coding':
							dftmp = dftmp[dftmp.Feature != 'exon']
						else:
							if len(dftmp[ dftmp.Feature.str.contains('utr') ]) > 0:
								# tr with bad annotation
								trRemove.append(tr)
						dftmp = dftmp.sort_values(by=['Start'])
						dftmp = dftmp.reset_index(drop=True)
						location = mapG4OnTr(dftmp, pG4r)
						pG4rtmp = pG4r
						pG4rtmp['Location'] = location
						pG4rtmp['Biotype'] = dftmp.Biotype[0]
						pG4rtmp.id[0] = tr
						dfpG4Annotation = dfpG4Annotation.append(pG4rtmp)
		else:
			pG4rtmp = dfpG4.iloc[pG4]
			id = pG4rtmp.id
			dftmpIntron = dfIntron[ dfIntron.Id == id ]
			dftmpIntron = dftmpIntron.reset_index(drop=True)
			pG4rtmp = mapG4onJunction(pG4rtmp, dftmpIntron, dfTr)
			dfpG4Annotation = dfpG4Annotation.append(pG4rtmp)
	dfpG4Annotation = dfpG4Annotation.reset_index(drop=True)
	print len(dfpG4Annotation)
	dfpG4Annotation = removeG4OnBadTr(dfpG4Annotation, trRemove)
	dfpG4Annotation = dfpG4Annotation.reset_index(drop=True)
	return dfpG4Annotation

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	sp = arg.specie
	filename = path+'/'+sp+'/'+sp+'.gtf'
	dfTr = Parser_gtf.importGTFdf(filename)
	dicoGene = Parser_gtf.importGTFGene(filename)
	main(dfTr, dicoGene, dfpG4)
