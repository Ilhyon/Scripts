#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

"""
``init.py`` **module description**:
This module has as input a fasta file of an entire chromosome but also
a gtf file from ensemble containing all inforation about genes and 
transcripts.
This script will parse the fasta file of a chromosome into a fasta file 
containing all genes of this chromosome.
For the assembly that were used, they are the one from pan-compara.
.. moduleauthor:: Anaìs Vannutelli, Jean-Pierre Perreault and Aida Ouangraoua
December 2018
Université de Sherbrooke Canada
Laboratoty CoBiUS and Jean-Pierre Perreault
"""

import re
import os
import math
import argparse
from pprint import pprint

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_Fasta')
	parser.add_argument ('-sp', '--specie', default = 'MM')
	return parser

def reverseSequence(Sequence) :
	"""
		Create the reverse complement of a DNA sequence.
	"""
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
			tmp = n # in some sequences there is many N or other letter
		reverse = reverse + tmp
	Sequence = reverse[::-1]
	return Sequence

def importFastaChromosome(sp) :
	"""
		Import the fasta of a chromosome, they are downloaded from 
		ensembl FTP for each specie. The assembly are those that are
		present in pan-compara.
		{chromosome : sequence}
	"""
	directory = "/home/anais/Documents/Data/Genomes/"+sp+"/Fasta/"
	dicoChromosome = {}
	listFile = os.listdir(directory)
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
	return dicoChromosome

def importGTF(sp, dicoChromosome) :
	"""
		Import information from the gtf file. The sequence is also added
		from the dictionary of chromosome.
		{idGene : {	Chromosome,
					Start,
					End,
					Strand,
					Sequence}}
	"""
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
							strand = "1"
						elif strand == "-":
							strand = "-1"
						if chrm in dicoChromosome :
							if int(startFeature) > 0 and endFeature > 0:
								"""
								-1 for coords because in python : start at 0 and not 1
								+1 for the sequence :
								1 2 3 4 5 6 7 8 9
								A T C G G G T A G
								exemple : 6-3 = 3 but there is 4 nucleotides
								so 6-3 +1 = 4
								"""
								geneSequence = dicoChromosome[chrm][(int(startFeature)-1):(int(endFeature)-1+1)]
							elif int(startFeature) < 0 :
								# genes overlaping the origin of 
								# replication in bacteria
								negSeq = dicoChromosome[chrm][int(startFeature):] # sequence before the origin of replication (for some bacteria)
								posSeq = dicoChromosome[chrm][0:int(endFeature)] # sequence after the origin of replication
								geneSequence = negSeq + posSeq
							dicoGene.update({idGene : {"Chromosome" : chrm, "Start" : startFeature, "End" : endFeature, "Strand" : strand, "Sequence" : geneSequence}})
	return dicoGene
	
def createFastaGene(sp, dicoGene):
	"""
		Create the output file, a fasta file with all gene of the specie.
	"""
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
			print "Ah"
		if strand == "-1" :
			Sequence = reverseSequence(Sequence)
		for i in range(0,int(nbLine)) :
			output.write(Sequence[cpt1:cpt2]+"\n") # to get the "real" fasta format
			cpt1 += 60
			cpt2 += 60
	output.close()
	
def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp=arg.specie # specie to analyse
	print "Fasta for "+sp
	dicoChromosome = importFastaChromosome(sp)
	dicoGene = importGTF(sp, dicoChromosome)
	createFastaGene(sp, dicoGene)
	print "\tDone"

main()
