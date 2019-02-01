#!/usr/bin/env python
# -*- coding: utf-8 -*-:


import argparse
from pprint import pprint

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_Fasta')
	parser.add_argument ('-s', '--specie', default = 'MM')
	return parser

def lengthG4() :
	directory = "/home/anais/Documents/Data/Mouse/results/MM_All_G4InTranscript.txt"
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
	directory = "/home/anais/Documents/Data/Mouse/All/MM_All_TranscriptType.txt"
	dicoTr = {}
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if l:
				words = l.split('|')
				idTr = words[0]
				coord = words[5].split(";")
				strand = words[3]
				if strand == "1":
					start = int(coord[0].split('-')[0])
					end = int(coord[-1].split('-')[1])
				else:
					end = int(coord[0].split('-')[0])
					start = int(coord[-1].split('-')[1])
				dicoTr[idTr] = {"Start" : start, "End" : end}
	return dicoTr
	
def lengthTranscriptUsed(dicoTr) :
	directory = "/home/anais/Documents/Data/Mouse/All/transcript_used.txt"
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















