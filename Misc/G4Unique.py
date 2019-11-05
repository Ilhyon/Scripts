#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import csv
import argparse
from pprint import pprint

########################################################################
def CreateListeG4Unique(filename):
	G4Unique = [] # list which will contain the unique G4s
	G4UniqueID = [] # just the position to know if it is alrady found or not
	with open(filename) as f: # file opening
		content = f.read()
		#~ print(content)
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			col = l.split('\t')
			InfoG4 = col[0].split('|')
			if len(InfoG4)>1: # if the line is not empty
				if InfoG4[1] not in G4UniqueID: # if it's the first time we encounter this G4s
				#we add it to the list
					#~ print(InfoG4[1])
					G4UniqueID.append(InfoG4[1])
					G4Unique.append(l)
		#~ print(len(G4Unique))
	return G4Unique
########################################################################
def ExtractionG4InTranscript(directory, specie, chromosome, G4InTranscript):
	output= open(directory+"/"+specie+"_"+chromosome+"_G4InTranscript_Unique.txt","w") ## file opening
	output.write("InfoG4ByTranscript\tcGcC\tG4Hunter\tsequenceG4\tG4NN\tlocalisation\ttranscriptBiotype\n")
	for i in G4InTranscript:
		output.write(i+"\n")

########################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Unique')
	parser.add_argument ('-p', '--path', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/Human/All/')
	parser.add_argument ('-CHR', '--CHROMOSOME', default = 'All')
	parser.add_argument ('-specie', '--specie', default = 'HS')
	return parser
########################################################################
def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path	# directory which contain all the directory chromosome
	CHROMOSOME=arg.CHROMOSOME	# chromosome to analyze
	specie=arg.specie	# specie to analyse
	
	fileG4 = path+specie+"_"+CHROMOSOME+"_G4InTranscript.txt"
	UniqueG4 = CreateListeG4Unique(fileG4)
	#~ print(UniqueG4)
	ExtractionG4InTranscript(path, specie, CHROMOSOME, UniqueG4)
	
########################################################################

main()
