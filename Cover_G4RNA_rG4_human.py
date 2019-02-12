#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import argparse
import re
from pprint import pprint

# This script allow to compute the number of pG4 that are present
# in the data base G4RNA. To do that we import the G4 from G4RNA. Only
# the folding G4 and the WT are imported. Then the pG4 are browsed. A 
# filter is done using the chromosome. Only the G4 (from G4RNA) that
# are in the same chromosome than the pG4 will be checked. Then the
# coverage is compute if their is an overlap between the position of 
# the pG4 and the G4, but also if the sequence match.
# The script will print the number of G4 that are imported, and then
# the number of pG4 that correspond to a G4.
# !! If many pG4 match a G4, it will be acounted only one time !!

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_Fasta')
	parser.add_argument ('-s', '--specie', default = 'MM')
	return parser

def importG4RNA():
	"""
		Read a file that contain all G4 from G4RNA screener (file from 31/01/2019)
		The db maybe update later so the results could change.
		Return a dictionary of G4 folded and WT : {chr : {start-end|length : {"Start","End", "Sequence"}}
	"""
	directory = "/home/anais/Documents/Data/G4RNA/G4rna_G4.txt"
	dicoG4RNA = {}
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		cpt = 0
		for l in lines:
			if not l.startswith('Gene') and l:
				words = l.split('\t')
				if words[2].split("r")[0] :
					chrm = words[2].split("r")[1]
				else :
					chrm = ""
				geneID = words[0]
				G4type = words[1]
				length = words[5]
				sequence = words[6]
				folding = words[7]
				start = words[3]
				end = words[4]
				if folding == "1" and geneID != "Artificial" and re.search("WT", G4type):
					if chrm not in dicoG4RNA :
						dicoG4RNA[chrm] = {str(start)+"-"+str(end) : {"Sequence" : sequence, "Start" : start, "End":end}}
					else :
						dicoG4RNA[chrm].update({str(start)+"-"+str(end) : {"Sequence" : sequence, "Start" : start, "End":end}})
	for chrm in dicoG4RNA :
		cpt += len(dicoG4RNA[chrm])
	print "Nombre de G4 dans G4RNA : "+str(cpt)
	return dicoG4RNA

def coverage(dicoG4RNA):
	"""
		Input : the dictionary of G4 
		Return : a list of G4 that are find in G4RNA and in pG4, I used 
		it to find the difference between the prediction using all score
		or only G4nn.
		It is possible to find if the G4 is a canonical G4 by uncomment 
		some lines, the sequences will be printed.
	"""
	# ~ directory = '/home/anais/Documents/Data/Human/G4Conserve-master/results/perChromosome/HS_All_G4InTranscript.txt' # only G4nn human
	directory = '/home/anais/Documents/Data/Human/All/HS_All_G4InTranscript.txt' # all score human
	e1 = r"(?i)(?P<run>g{3,})(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)" # regular expression to find canonical G4
	result = []
	
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if not l.startswith('Info') and l:
				words = l.split('\t')
				sequence = words[3]
				header = words[0]
				transcriptID = header.split("|")[0]
				idG4 = header.split("|")[1]
				chrm = idG4.split(":")[0]
				start = idG4.split(":")[1].split("-")[0]
				end = idG4.split(":")[1].split("-")[1]
				if chrm in dicoG4RNA :
					for G4 in dicoG4RNA[chrm]:
						coord = G4.split("|")[0]
						startG4 = coord.split("-")[0]
						endG4 = coord.split("-")[1]
						G4seq = dicoG4RNA[chrm][G4]["Sequence"]
						if G4+"\t"+G4seq not in result :
							if (startG4 > start and startG4 < end) or (endG4 > start and endG4 < end):
								result.append(G4+"\t"+G4seq)
								# ~ if(re.search(e1,l)):
									# ~ print re.findall(e1,sequence)
						if G4+"\t"+G4seq not in result :
							if (start > startG4 and start < endG4) or (end > startG4 and end < endG4):
								result.append(G4+"\t"+G4seq)
								# ~ if(re.search(e1,l)):
									# ~ print re.findall(e1,sequence)
						if G4+"\t"+G4seq not in result :
							if (re.search(G4seq, sequence) or re.search(sequence, G4seq)):
								# ~ if(re.search(e1,l)):
									# ~ print re.findall(e1,sequence)
								result.append(G4+"\t"+G4seq)
	print "Nombre de G4 de G4RNA trouvés dans les pG4 : "+str(len(result))
	return result
	
def writepG4(listpG4):
	output= open("/home/anais/Documents/Data/G4RNA/G4_all_score","w")
	output.write("\n".join(listpG4))


def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp=arg.specie	# specie to analyse
	dicoG4RNA = importG4RNA()
	listpG4 = coverage(dicoG4RNA)
	writepG4(listpG4)

main()
