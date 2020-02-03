#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
``init.py`` **module description**:
This module has as input the output of G4RNAScreener and a dictionnary
with all the parameters that were used for G4RNA Screener.
This script filter the output of G4RNA Screener to only keep windows
over thresholds. Then all those windows are merged if they are
overlapping. Overlapping windows are pG4.
.. moduleauthor:: Anaìs Vannutelli, Michel Garant, Sarah Bellamiti, Jean-Pierre Perreault and Aida Ouangraoua
March 2019
Université de Sherbrooke Canada
Laboratoty CoBiUS and Jean-Pierre Perreault
"""

import os
import getpG4
import argparse
import Parser_gtf
import G4Annotation
import pandas as pd
from pprint import pprint
import recurrentFunction as rF

def importIntron(filename):
	try:
		df = pd.read_csv(filename, sep='\t')
	except:
		print("This file couldn't be converted in data frame : " + filename)
	else:
		# dataFrame with all windows from G4RNA Screener
		df.columns = ['Transcript', 'Chromosome','Start','End', 'Strand']
		df['Id'] = df['Start'].map(str) +':'+ df['End'].map(str)
		return df

def mergeWindow(path, dicoParam, option):
	if option == 'Annotation' or option == 'venn':
		directory = path + '/CSVFile'
	elif option == 'Random':
		directory = path + '/randomCSV'
	dfpG4 = pd.DataFrame()
	dicoVenn = {}
	# directory containing data for a specie
	for path, dirs, files in os.walk(directory):
		# for each element of the directory to passed
		for filename in files: # for each files
			inputfile = directory + '/' + filename
			if option == 'Random' or option == 'Annotation':
				if ('gene_unspliced' in filename):
					# windows in genes
					dfpG4 = dfpG4.append(getpG4.main(inputfile,
							dicoParam, "Gene", option))
					dfpG4 = dfpG4.reset_index(drop=True)
				elif ('transcript_unspliced' in filename):
					dfpG4 = dfpG4.append(getpG4.main(inputfile,
							dicoParam,"Junction", option))
					dfpG4 = dfpG4.reset_index(drop=True)
				elif ('Random' in filename):
					dfpG4 = dfpG4.append(getpG4.main(inputfile,
							dicoParam,"Gene", option))
					dfpG4 = dfpG4.reset_index(drop=True)
			elif option == 'venn':
				if ('gene_unspliced' in filename):
					dicoTmp = getpG4.mainControl(inputfile, dicoParam, 'Gene', option)
				else:
					dicoTmp = {}
				for score in dicoTmp:
					if score in dicoVenn:
						dicoVenn[score].extend(dicoTmp[score])
					else:
						dicoVenn[score] = dicoTmp[score]
	if len(dfpG4) > 0:
		dfpG4 = dfpG4.drop_duplicates(subset=None, keep='first', inplace=False)
		dfpG4 = dfpG4.reset_index(drop=True)
	if option == 'Random' or option == 'Annotation':
		return dfpG4
	elif option == 'venn':
		return dicoVenn

def main(dicoParam, path, dicoGene, dicoTr, dicoUTR, dfIntron, option):
	if option == 'Annotation':
		output1 = path + '/pG4.txt'
		output2 = path + '/pG4Mature.txt'
		output0 = path + '/G4.txt'
		dfpG4 = mergeWindow(path, dicoParam, option)
		dfpG4.to_csv(path_or_buf=output0, header=True, index=None, sep='\t')
		print ('\t - Number of pG4 in genes : ' + str(len(dfpG4.index)))
		# print(len(dfpG4[dfpG4.Feature == 'Gene']))
		pG4Anno, pG4MatureAnno = G4Annotation.main(dicoTr, dicoGene, dicoUTR, dfpG4, dfIntron)
		pG4Anno = pG4Anno.drop_duplicates(subset=None, keep='first', inplace=False)
		pG4Anno.to_csv(path_or_buf=output1, header=True, index=None, sep='\t')
		pG4MatureAnno.to_csv(path_or_buf=output2, header=True, index=None, sep='\t')
		print ('\t - Number of pG4 in preRNA : ' + str(len(pG4Anno.index)))
		print ('\t - Number of pG4 in matureRNA : ' + str(len(pG4MatureAnno.index)))
	elif option == 'Random':
		dfpG4 = mergeWindow(path, dicoParam, option)
		print ('\t'+str(dfpG4.shape))
		print(dfpG4)


def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	parser.add_argument ('-sp', '--specie', default = \
		'yersinia_pestis_biovar_microtus_str_91001')
	parser.add_argument ('-o', '--option', default = 'Annotation')
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
	sp = arg.specie
	option = arg.option
	ini = rF.setUpperLetter(sp)
	path = arg.path + sp
	print("Specie : " + sp)
	dicoParam = rF.createDicoParam(arg)
	# dfTr = Parser_gtf.importGTFdf(path +'/'+ sp +'.gtf')
	dicoTr = Parser_gtf.importTranscriptFromGTF(path +'/'+ sp +'.gtf')
	dicoUTR = Parser_gtf.importUTRFromGTF(path +'/'+ sp +'.gtf')
	dicoGene = Parser_gtf.extractDicoGeneTr(path +'/'+ sp +'.gtf')
	dfIntron = importIntron(path +'/'+ ini +'_intron.txt')
	main(dicoParam, path, dicoGene, dicoTr, dicoUTR, dfIntron, option)
	print("\tDone")
