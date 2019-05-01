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

def main(dicoTr, dicoGene, dfpG4):
	cpt = 0
	for tr in dicoTr:
		trStart = dicoTr[tr]["Start"]
		trEnd = dicoTr[tr]["End"]
		gene = dicoTr[tr]["Gene"]
		dfG4Ingene = dfpG4[dfpG4.Gene == gene]
		if not dfG4Ingene.empty :
			cpt += 1
	print cpt

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	sp = arg.specie
	filename = path+'/'+sp+'/'+sp+'.gtf'
	dicoTr = Parser_gtf.importGTF(filename)
	dicoGene = Parser_gtf.importGTFGene(filename)
	main(dicoTr, dicoGene, dfpG4)
