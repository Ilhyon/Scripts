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
		Return a dictionary of G4 folded and WT : {chr : {start-end : {"Start","End", "Sequence"}}
	"""
	directory = "/home/anais/Documents/Data/G4RNA/G4rna_G4.txt"
	dicoG4RNA = {}
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		cpt = 0
		list_G4_to_Delete = []
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
				sequence = words[7]
				folding = words[8]
				start = words[3]
				end = words[4]
				strand = words[6]
				if strand == "+":
					strand = "1"
				else :
					strand = "-1"
				if folding == "1" and geneID != "Artificial" and re.search("WT", G4type):
					if chrm not in dicoG4RNA :
						dicoG4RNA[chrm] = {str(start)+"-"+str(end) : {"Sequence" : sequence,
																	"Start" : start,
																	"End" : end,
																	"Length" : length,
																	"Strand" : strand}}
					else :
						for G4 in dicoG4RNA[chrm] :
							start_dico = G4.split("-")[0]
							end_dico = G4.split("-")[1]
							if (start == start_dico or end == end_dico) and length > dicoG4RNA[chrm][G4]["Length"]:
								list_G4_to_Delete.append([chrm, G4])
						dicoG4RNA[chrm].update({str(start)+"-"+str(end) : {"Sequence" : sequence,
																			"Start" : start,
																			"End":end,
																			"Length" : length,
																			"Strand" : strand}})
	for G4_to_Delete in list_G4_to_Delete :
		chrm = G4_to_Delete[0]
		coordG4 = G4_to_Delete[1]
		if coordG4 in dicoG4RNA[chrm] :
			del dicoG4RNA[chrm][coordG4]
	for chrm in dicoG4RNA :
		cpt += len(dicoG4RNA[chrm])
	# ~ pprint(dicoG4RNA)
	print "Number of G4 in G4RNA : "+str(cpt)
	return dicoG4RNA

def coveragepG4(dicoTarget, Dtype, directory, bd):
	"""
		Input : the dictionary of G4
	"""
	resultrG4 = []
	resultpG4 = []
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if not l.startswith('Info') and l:
				words = l.split('\t')
				sequence = words[3]
				header = words[0]
				transcriptID = header.split("|")[0]
				pG4 = header.split("|")[1]
				strandpG4 = header.split("|")[2]
				chrm = pG4.split(":")[0]
				startpG4 = pG4.split(":")[1].split("-")[0]
				endpG4 = pG4.split(":")[1].split("-")[1]
				biotype = words[6]
				location = words[5]
				if location != "Intron" and biotype == "protein_coding" :
					if chrm in dicoTarget :
						for G4 in dicoTarget[chrm]:
							if strandpG4 == dicoTarget[chrm][G4]["Strand"]:
								coord = G4.split("|")[0]
								startG4 = coord.split("-")[0]
								endG4 = coord.split("-")[1]
								G4seq = dicoTarget[chrm][G4]["Sequence"]
								if G4+"\t"+G4seq not in resultrG4 :
									if (startG4 > startpG4 and startG4 < endpG4) or (endG4 > startpG4 and endG4 < endpG4):
										resultrG4.append(G4+"\t"+G4seq)
										resultpG4.append(pG4+"\t"+sequence)
								if G4+"\t"+G4seq not in resultrG4 :
									if (startpG4 > startG4 and startpG4 < endG4) or (endpG4 > startG4 and endpG4 < endG4):
										resultrG4.append(G4+"\t"+G4seq)
										resultpG4.append(pG4+"\t"+sequence)
								if G4seq not in resultrG4 :
									if (re.search(G4seq, sequence) or re.search(sequence, G4seq)):
										resultrG4.append(G4seq)
										resultpG4.append(pG4+"\t"+sequence)
	print "Number of pG4r ("+Dtype+") corresponding to G4 from "+bd+" : "+str(len(resultrG4))
	print "Number of pG4r ("+Dtype+") corresponding to G4 from "+bd+" : "+str(len(list(set(resultpG4))))
	return resultrG4

def coveragerG4seq(dicoG4RNA, Dtype, directory):
	"""
		Input : the dictionary of G4 
		Return : a list of G4 that are find in G4RNA and in pG4, I used 
		it to find the difference between the prediction using all score
		or only G4nn.
	"""
	result = []
	DicorG4seq = {}
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if not l.startswith('#') and l:
				words = l.split('\t')
				sequence = words[9]
				transcriptID = words[10]
				chrm = words[0].split("r")[1]
				startrG4seq = words[1]
				endrG4seq = words[2]
				strandrG4 = words[5]
				if strandrG4 == "+":
					strandrG4 = "1"
				else :
					strandrG4 = "-1"
				if chrm in DicorG4seq :
					DicorG4seq[chrm].update({str(startrG4seq)+"-"+str(endrG4seq) : {"Sequence" : sequence,
																					"Start" : startrG4seq,
																					"End":endrG4seq,
																					"Strand":strandrG4}})
				else :
					DicorG4seq.update({chrm : {str(startrG4seq)+"-"+str(endrG4seq) : {"Sequence" : sequence,
																					"Start" : startrG4seq,
																					"End":endrG4seq,
																					"Strand":strandrG4}}})
				if chrm in dicoG4RNA :
					for G4 in dicoG4RNA[chrm]:
						coord = G4.split("|")[0]
						startG4 = coord.split("-")[0]
						endG4 = coord.split("-")[1]
						G4seq = dicoG4RNA[chrm][G4]["Sequence"]
						if G4+"\t"+G4seq not in result :
							if (startG4 > startrG4seq and startG4 < endrG4seq) or (endG4 > startrG4seq and endG4 < endrG4seq):
								result.append(G4+"\t"+G4seq)
						if G4+"\t"+G4seq not in result :
							if (startrG4seq > startG4 and startrG4seq < endG4) or (endrG4seq > startG4 and endrG4seq < endG4):
								result.append(G4+"\t"+G4seq)
						if G4+"\t"+G4seq not in result :
							if (re.search(G4seq, sequence) or re.search(sequence, G4seq)):
								result.append(G4+"\t"+G4seq)
	print "Number of "+Dtype+" corresponding to G4 from G4RNA : "+str(len(result))
	return result, DicorG4seq

def findCanonicalG4(liste) :
	"""
		Test for each G4 from G4RNA that were find in the predictions 
		if it is a canonical G4 or not. 
		The function doesn't return anything but print the number of 
		canonical and non canonical G4.
	"""
	nbCanonical = 0
	e1 = r"(?i)(?P<run>g{3,})(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)" # regular expression to find canonical G4
	for G4 in liste :
		coord = G4.split("\t")[0]
		sequence = G4.split("\t")[1]
		if(re.search(e1, sequence)):
			nbCanonical += 1
	print "\tAmong those there is "+str(nbCanonical)+" canonicals G4 and "+str(len(liste)-nbCanonical)+" non canonical G4"

def findCanonicalG4FromDico(dico) :
	"""
		Test for each G4 from G4RNA if it is a canonical G4 or not. 
		The function doesn't return anything but print the number of 
		canonical and non canonical G4.
	"""
	nbCanonical = 0
	nbTot = 0
	e1 = r"(?i)(?P<run>g{3,})(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)" # regular expression to find canonical G4
	for chrm in dico :
		nbTot += len(dico[chrm])
		for G4 in dico[chrm] :
			sequence = dico[chrm][G4]["Sequence"]
			if(re.search(e1, sequence)):
				nbCanonical += 1
	print "\tAmong those there is "+str(nbCanonical)+" canonicals G4 and "+str(nbTot-nbCanonical)+" non canonical G4"

def writepG4(listpG4):
	output= open("/home/anais/Documents/Data/G4RNA/G4_all_score","w")
	output.write("\n".join(listpG4))


def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp=arg.specie	# specie to analyse
	directoryG4nn = '/home/anais/Documents/Data/Human/G4Conserve-master/results/perChromosome/HS_All_G4InTranscript.txt' # only G4nn human
	directoryAll = '/home/anais/Documents/Data/Human/All/HS_All_G4InTranscript.txt' # all score human
	directoryrG4seq_K = '/home/anais/Documents/Publis/rG4seq/cdt_K.csv'
	directoryrG4seq_PDS_K = '/home/anais/Documents/Publis/rG4seq/cdt_PDS-K.csv'
	dicoG4RNA = importG4RNA()
	# ~ findCanonicalG4FromDico(dicoG4RNA)
	# ~ listpG4G4nn = coveragepG4(dicoG4RNA, "G4nn", directoryG4nn, "G4RNA")
	# ~ findCanonicalG4(listpG4G4nn)
	# ~ listpG4All = coveragepG4(dicoG4RNA, "all score", directoryAll, "G4RNA")
	# ~ findCanonicalG4(listpG4All)
	listpG4rG4seq_K, dicoK = coveragerG4seq(dicoG4RNA, "rG4seq K+", directoryrG4seq_K)
	# ~ findCanonicalG4(listpG4rG4seq_K)
	listrG4K = coveragepG4(dicoK, "G4all", directoryAll, "Kwok K+")
	# ~ findCanonicalG4(listrG4K)
	# ~ listpG4rG4seq_K, dicoG4seq_PDS_K = coveragerG4seq(dicoG4RNA, "rG4seq K+", directoryrG4seq_PDS_K)
	# ~ listrG4K = coveragepG4(dicoG4seq_PDS_K, "G4all", directoryAll, "Kwok K+")
	del dicoK
	
main()
