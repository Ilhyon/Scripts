#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import re
import argparse

########################################################################
def parserBold(filename):
	fileparse = []
	e1 = r"(?i)(?P<run>g{3,})(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)"
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			#~ print "---------------------------------"
			#~ print l
			if(re.search(e1,l)):
				print re.findall(e1,l)
########################################################################
def ExtractionG4InTranscript(directory, specie, chromosome, G4InTranscript):
	output= open(directory+"/"+specie+"_chr"+chromosome+"_G4InTranscript_Unique.txt","w") ## file opening
	output.write("InfoG4ByTranscript\tcGcC\tG4Hunter\tsequenceG4\tG4NN\tlocalisation\ttranscriptBiotype\n")
	for i in G4InTranscript:
		output.write(i+"\n")
########################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Unique')
	parser.add_argument ('-p', '--path', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/Human/All/non_stop_decay.txt')
	parser.add_argument ('-sf', '--subfamily', default = 'non_stop_decay')
	parser.add_argument ('-s', '--specie', default = 'MM')
	return parser
########################################################################
def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path	# directory which contain all the directory chromosome
	specie=arg.specie	# specie to analyse
	subfamily=arg.subfamily 
	
	fileSubf = path
	parserBold(fileSubf)	
	
########################################################################

main()
