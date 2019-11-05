#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import re
import argparse

########################################################################
def parserBold(filename):
	e1 = r"(?i)(?P<run>g{3,})(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)"
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if(re.search(e1,l)):
				print l
			# ~ if not l.startswith('Gene') and l:
				# ~ words = l.split('\t')
				# ~ if words[2].split("r")[0] :
					# ~ chrm = words[2].split("r")[1]
				# ~ else :
					# ~ chrm = ""
				# ~ geneID = words[0]
				# ~ sequence = words[6]
				# ~ G4type= words[1]
				# ~ folding = words[7]
				# ~ if folding == "1" and geneID != "Artificial" and re.search("WT", G4type):
					# ~ if(re.search(e1,l)):
						# ~ print re.findall(e1,sequence)
########################################################################
def ExtractionG4InTranscript(directory, specie, chromosome, G4InTranscript):
	output= open(directory+"/"+specie+"_chr"+chromosome+"_G4InTranscript_Unique.txt","w") ## file opening
	output.write("InfoG4ByTranscript\tcGcC\tG4Hunter\tsequenceG4\tG4NN\tlocalisation\ttranscriptBiotype\n")
	for i in G4InTranscript:
		output.write(i+"\n")
########################################################################
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Unique')
	parser.add_argument ('-p', '--path', default = '/home/anais/Documents/Data/G4RNA/G4_G4nn')
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
