#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import argparse
import re
from pprint import pprint

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_Fasta')
	parser.add_argument ('-s', '--specie', default = 'MM')
	return parser

def importG4RNA():
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
					cpt +=1
					if geneID not in dicoG4RNA :
						dicoG4RNA[geneID] = {str(start)+"-"+str(end)+"|"+str(length) : {"Sequence" : sequence, "Chromosome" : chrm, "Start" : start, "End":end}}
					else :
						dicoG4RNA[geneID].update({str(start)+"-"+str(end)+"|"+str(length) : {"Sequence" : sequence, "Chromosome" : chrm, "Start" : start, "End":end}})
	print "Nombre de G4 dans G4RNA : "+str(cpt)
	return dicoG4RNA

def importID():
	directory = "/home/anais/Documents/Data/G4RNA/EnsemblID_to_HGNCSymbols.txt"
	dicoID = {}
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if not l.startswith('Gene') and l:
				words = l.split('\t')
				geneID = words[0]
				transcriptID = words[1]
				hgncID = words[2]
				dicoID[transcriptID] = hgncID
	return dicoID

def Coverage(dicoG4RNA, dicoID):
	# ~ directory = '/home/anais/Documents/Data/Human/G4Conserve-master/results/perChromosome/HS_All_G4InTranscript.txt' # only G4nn human
	directory = '/home/anais/Documents/Data/Human/All/HS_All_G4InTranscript.txt' # all score human
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
				if transcriptID in dicoID : # if the id ensembl of the transcript exist in the dictionary
					if dicoID[transcriptID] in dicoG4RNA : # if the hgnc id have a G4 in G4rna
						hgncID = dicoID[transcriptID]
						for G4 in dicoG4RNA[hgncID]:
							coord = G4.split("|")[0]
							startG4 = coord.split("-")[0]
							endG4 = coord.split("-")[1]
							chrmG4 = dicoG4RNA[hgncID][G4]["Chromosome"]
							G4seq = dicoG4RNA[hgncID][G4]["Sequence"]
							if G4 not in result :
								if chrm == chrmG4 :
									if G4 not in result :
										if (startG4 > start and startG4 < end) or (endG4 > start and endG4 < end):
											result.append(G4)
									if G4 not in result :
										if (start > startG4 and start < endG4) or (end > startG4 and end < endG4):
											result.append(G4)
							if G4 not in result :
								if (re.search(G4seq, sequence) or re.search(sequence, G4seq)):
									result.append(G4)
	return result


def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp=arg.specie	# specie to analyse
	# ~ directory = "/home/anais/Documents/Data/Genomes/"+sp+"/Fasta/"
	dicoG4RNA = importG4RNA()
	# ~ pprint(dicoG4RNA)
	dicoID = importID()
	print "Nombre de G4 de G4RNA trouvés dans les pG4 : "+str(len(Coverage(dicoG4RNA,dicoID)))
	

main()
