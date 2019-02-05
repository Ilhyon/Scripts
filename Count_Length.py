#!/usr/bin/env python
# -*- coding: utf-8 -*-:


import argparse
from pprint import pprint

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_Fasta')
	parser.add_argument ('-s', '--specie', default = 'MM')
	return parser

def lengthG4() :
	directory = "/home/anais/Documents/Data/Mouse/mouseEssai/out_all_score/MM_All_G4InTranscript.txt"
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
	return totLengthG4

def importInfoTr():
	directory = "/home/anais/Documents/Data/Mouse/All/MM_All_Transcript_location.txt"
	dicoTr = {}
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if l:
				words = l.split('\t')
				idTr = words[0]
				strand = words[1]
				# ~ exonList = words[5]
				# ~ exonList = exonList.split(";")
				# ~ exon_total = []
				# ~ for i in exonList :
					# ~ if strand == 1 :
						# ~ exon_total.append([int(i.split("-")[0]),int(i.split("-")[1])])
					# ~ else :
						# ~ exon_total.append([int(i.split("-")[1]),int(i.split("-")[0])])
				# ~ exonSorted = []
				# ~ for exon in exon_total:
					# ~ if not exonSorted :
						# ~ exonSorted = [exon]
					# ~ elif exon[0] < exonSorted[0][0]: # the new exon start is the first exon
						# ~ exonSorted.insert(0, exon) # then we change the first exon
					# ~ elif exon[0] > exonSorted[-1][0]: # the new exon start is the last exon
						# ~ exonSorted.insert(-1, exon)
					# ~ else:
						# ~ for i in exon_total[1::-2]:
							# ~ if exon[0] < i[0]:
								# ~ exonSorted.insert(exon_total.index(i),exon)
							# ~ else:
								# ~ exonSorted.insert(exon_total.index(i)+1,exon)
				start = int(words[2])
				end = int(words[3])
				dicoTr[idTr] = {"Start" : start , "End" : end}
	return dicoTr
	
def lengthTranscriptUsed(dicoTr) :
	directory = "/home/anais/Documents/Data/Mouse/All/transcript_used_MM.txt"
	totLengthTr = 0
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if l:
				words = l.split('\t')
				idTr = words[0]
				if idTr in dicoTr:
					totLengthTr += dicoTr[idTr]["End"] - dicoTr[idTr]["Start"] + 1
	return totLengthTr
	
def lengthTranscriptome(dicoTr) :
	totLengthTr = 0
	for tr in dicoTr:
		totLengthTr += dicoTr[tr]["End"] - dicoTr[tr]["Start"] + 1
	return totLengthTr

def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp=arg.specie	# specie to analyse
	# ~ directory = "/home/anais/Documents/Data/Genomes/"+sp+"/Fasta/"
	print "Length of rG4 : "+str(lengthG4())
	dicoTr = importInfoTr()
	print "Lenght of transcript used : "+str(lengthTranscriptUsed(dicoTr))
	print "Length of the entire transcriptome : "+str(lengthTranscriptome(dicoTr))

main()















