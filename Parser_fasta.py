#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import os
import math
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
	directory = "/home/anais/Documents/Data/Genomes/"+sp+"/Fasta/"
	dicoChromosome = {}
	listFile = os.listdir(directory)
	print "Fasta for "+sp
	for filename in listFile :
		with open(directory+filename) as f: # file opening
			content = f.read()
			l = content.split('\n')
			if l[0].startswith('>'):
				header = l[0]
				chrm = header.split(' ')[0][1:]
				sequence = "".join(l[1:])
				dicoChromosome.update({chrm : sequence})
	print "dico fasta done"
	filename = "/home/anais/Documents/Data/Genomes/"+sp+"/"+sp+".gtf"
	exists = os.path.isfile(filename)
	if exists :	
		dicoGene = {}
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					words=l.split('\t')
					idGene = words[8].split(';')[0].split('"')[1]
					feature = words[2]
					if feature == "gene" :
						chrm = words[0]
						startFeature = words[3]
						endFeature = words[4]
						strand = words[6]
						if strand == "+":
							strand = 1
						elif strand == "-":
							strand = -1
						if chrm in dicoChromosome :
							geneSequence = dicoChromosome[chrm][(int(startFeature)-1):(int(endFeature)-1+1)]
							dicoGene.update({idGene : {"Chromosome" : chrm, "Start" : startFeature, "End" : endFeature, "Strand" : str(strand), "Sequence" : geneSequence}})
	
	words = sp.split("_")
	letters = [word[0] for word in words]
	ini = "".join(letters)
	ini = ini.upper()
	output = open("/home/anais/Documents/Data/Genomes/"+sp+"/"+ini+"_gene_unspliced.txt","w") # file opening for reading
	for gene in dicoGene :
		chromosome = dicoGene[gene]["Chromosome"]
		start = dicoGene[gene]["Start"]
		end = dicoGene[gene]["End"]
		Sequence = dicoGene[gene]["Sequence"]
		strand = str(dicoGene[gene]["Strand"])
		output.write(">"+gene+"|"+chromosome+"|"+strand+"|"+start+"|"+end+"\n")
		nbLine = math.ceil(float(int(end)-int(start))/60)
		cpt1 = 0
		cpt2 = 60
		if nbLine < 0 :
			print AH
		if strand == "-1" :
			reverse = ""
			for n in Sequence:
				if n == "A" :
					tmp = "T"
				elif n == "T" :
					tmp = "A"
				elif n == "G" :
					tmp = "C"
				elif n == "C" :
					tmp = "G"
				else :
					tmp = n
				reverse = reverse + tmp
			Sequence = reverse[::-1]
		for i in range(0,int(nbLine)) :
			output.write(Sequence[cpt1:cpt2]+"\n")
			cpt1 += 60
			cpt2 += 60
	output.close()
	print "Done"
	
# ~ pprint(dicoGene)

main()
