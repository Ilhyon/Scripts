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
	location = ''
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

def main(dfTr, dicoGene, dfpG4):
	dfpG4Annotation = pd.DataFrame()
	for pG4 in range(0,len(dfpG4)):
		geneId = dfpG4.id.iloc[pG4]
		if geneId in dicoGene:
			if 'Transcript' in dicoGene[geneId]:
				for tr in dicoGene[geneId]['Transcript']:
					pG4r = dfpG4.iloc[pG4:pG4+1]
					pG4r = pG4r.reset_index(drop=True)
					if pG4r.Feature[0] == 'Gene':
						dftmp =  dfTr[ dfTr.Transcript == tr ].dropna()
						if dftmp.Type.iloc[0] == 'Coding':
							dftmp = dftmp[dftmp.Feature != 'exon']
						dftmp = dftmp.sort_values(by=['Start'])
						dftmp = dftmp.reset_index(drop=True)
						location = mapG4OnTr(dftmp, pG4r)
						pG4rtmp = pG4r
						pG4rtmp['Location'] = location
						pG4rtmp['Biotype'] = dftmp.Biotype[0]
						pG4rtmp['id'][0] = tr
						dfpG4Annotation = dfpG4Annotation.append(pG4rtmp)
					else:
						pG4rtmp = pG4r
						pG4rtmp['Location'] = 'junction'
						pG4rtmp['Biotype'] = dftmp.Biotype[0]
						dfpG4Annotation = dfpG4Annotation.append(pG4rtmp)
	dfpG4Annotation = dfpG4Annotation.reset_index(drop=True)
	print dfpG4Annotation

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	sp = arg.specie
	filename = path+'/'+sp+'/'+sp+'.gtf'
	dfTr = Parser_gtf.importGTFdf(filename)
	dicoGene = Parser_gtf.importGTFGene(filename)
	main(dfTr, dicoGene, dfpG4)
