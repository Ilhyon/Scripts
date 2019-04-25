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

def createListCodingProtein():
	codingProtein=['IG_C_gene',
					'IG_D_gene',
					'IG_J_gene',
					'IG_LV_gene',
					'IG_M_gene',
					'IG_V_gene',
					'IG_Z_gene',
					'nonsense_mediated_decay',
					'nontranslating_CDS',
					'non_stop_decay',
					'protein_coding',
					'TR_C_gene',
					'TR_D_gene',
					'TR_gene',
					'TR_J_gene',
					'TR_V_gene']
	return codingProtein

def main(dicoParam, path):
	codingProtein = createListCodingProtein()
	dfpG4 = pd.DataFrame()
	G4DetectedInJunction = {}
	directory = path + '/CSVFile'
	# directory containing data for a specie
	for path, dirs, files in os.walk(directory):
		# for each element of the directory to passed
		for filename in files: # for each files
			inputfile = directory + '/' + filename
			if ('gene_unspliced' in filename):
				# windows in genes
				dfpG4 = dfpG4.append(getpG4.main(inputfile,
						dicoParam, "Gene"))
				dfpG4 = dfpG4.reset_index(drop=True)
			# elif ('transcript_unspliced' in filename and '.csv' in filename): ## for G4 in junction CDS-CDS --> from splicing
			# 	dfpG4 = dfpG4.append(getpG4.main(inputfile,
			# 			dicoParam,"Junction"))
			# 	dfpG4 = dfpG4.reset_index(drop=True)
	dfpG4 = dfpG4.drop_duplicates(subset=None, keep='first', inplace=False)
	print '\t'+str(len(dfpG4))
	with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
		print(dfpG4)
	# G4Annotation.main(dicoTr, dicoGene, dfpG4)

def createDicoParam(arg):
	dicoParam = {"g4H" : float(arg.THRESHOLD_G4H),
				"cGcC" : float(arg.THRESHOLD_CGCC),
				"g4NN" : float(arg.THRESHOLD_G4NN),
				"junctionLength" : int(arg.EXTENSION),
				"windowLength" : int(arg.EXTENSION),
				"step" : int(arg.STEP)}
	return dicoParam

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	parser.add_argument ('-sp', '--specie', default = \
		'yersinia_pestis_biovar_microtus_str_91001')
	parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-E', '--EXTENSION', default = 100)
	parser.add_argument ('-W', '--WINDOW', default = 60)
	parser.add_argument ('-S', '--STEP', default = 10)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp = arg.specie
	path = arg.path + sp
	dicoParam = createDicoParam(arg)
	print "Specie : " + sp
	main(dicoParam, path)
	print "\tDone"


	# ~ listeG4InGeneEntire={}
	# ~ listeG4InGeneJunction={}
