#!/usr/bin/env python
# -*- coding: utf-8 -*-:


import argparse
from pprint import pprint

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_Fasta')
	parser.add_argument ('-s', '--specie', default = 'MM')
	return parser
	
def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp=arg.specie	# specie to analyse
	# ~ directory = "/home/anais/Documents/Data/Genomes/"+sp+"/Fasta/"
	directory = "/home/anais/Documents/Data/Human/G4Conserve-master/results/perChromosome/HS_All_G4InTranscript.txt"
	totLengthG4 = 0
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if not l.startswith('InfoG4ByTranscript') and l:
				words = l.split('\t')
				idG4 = words[0]
				idG4 = idG4.split('|')[1]
				startG4 = idG4.split(':')[1].split('-')[0]
				endG4 = idG4.split(':')[1].split('-')[1]
				lengthG4 = int(endG4) - int(startG4) + 1
				totLengthG4 += lengthG4
	print "Length of rG4 : "+totLengthG4

main()
