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
.. moduleauthor:: Anaìs Vannutelli, Aida Ouangraoua
December 2018
Université de Sherbrooke Canada
Laboratoty CoBiUS and Jean-Pierre Perreault
"""

import re
import os
import math
import argparse
import Parser_gtf as pgtf
from pprint import pprint

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_Fasta')
	parser.add_argument ('-sp', '--specie', default = 'MM')
	return parser

def reverseSequence(Sequence) :
	""" Reverse complement a DNA sequence.

	:param Sequence: DNA sequence that will be reversed.
	:type Sequence: string

	:returns: Sequence, the initial DNA sequence but reverse complemented.
	:rtype: string
	"""
	reverse = ""
	nucleotides = {'A' : 'T',
					'T' : 'A',
					'C' : 'G',
					'G' : 'C'}
	for n in Sequence:
		if n in nucleotides:
			tmp = nucleotides[n]
		else :
			tmp = n # in some sequences there is many N or other letter
		reverse += tmp
	Sequence = reverse[::-1]
	return Sequence

def importFastaChromosome(sp) :
	""" Import the fasta file contain the whole chromosome sequence.

	Fasta files are downloaded from ensembl FTP for each specie. The assembly
	are those that are present in pan-compara.

	:param sp: current specie.
	:type sp: string

	:returns: dicoChromosome, all chromosome and their sequences.
	:rtype: dictionary
	"""
	directory = "/home/anais/Documents/Data/Genomes/" + sp + "/Fasta/"
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

def getGeneSequence(sp, dicoChromosome) :
	""" Retrieves the sequence of all genes from a specie.

	Import gtf informations about genes, and then retrives the sequence of all
	genes. In some bacteria a gene can overlap the origine of replication. In
	that case the start of the gene is negative.

	:param sp: current specie.
	:type sp: string
	:param dicoChromosome: all chromosome and their sequences.
	:type dicoChromosome: dictionary

	:returns: dicoGene, all genes and their sequences, coordinates, etc.
	:rtype: dictionary
	"""
	filename = "/home/anais/Documents/Data/Genomes/" + sp + "/" + sp + ".gtf"
	exists = os.path.isfile(filename)
	if exists :
		dicoGene = pgtf.importGTFSequence(filename)
		for geneId in dicoGene:
			if dicoGene[geneId]['Chromosome'] in dicoChromosome :
				if (dicoGene[geneId]['geneStart'] > 0 and
					dicoGene[geneId]['geneEnd'] > 0):
					"""
					-1 for coords because in python : start at 0 and not 1
					+1 for the sequence :
					1 2 3 4 5 6 7 8 9
					A T C G G G T A G
					exemple : 6-3 = 3 but there is 4 nucleotides
					so 6-3 +1 = 4
					"""
					geneSequence = \
						dicoChromosome[dicoGene[geneId]['Chromosome']]\
						[ dicoGene[geneId]['geneStart']-1:\
						(dicoGene[geneId]['geneEnd']-1+1) ]
				elif dicoGene[geneId]['geneStart'] < 0 :
					# genes overlaping the origin of
					# replication in bacteria
					negSeq = dicoChromosome[chrm]\
						[int(dicoGene[geneId]['Chromosome']):]
						# sequence before the origin of replication
					posSeq = dicoChromosome[chrm]\
						[0:int(dicoGene[geneId]['geneEnd'])]
						# sequence after the origin of replication
					geneSequence = negSeq + posSeq
				dicoGene[geneId].update({"Sequence" : geneSequence})
	return dicoGene

def createFastaGene(sp, dicoGene):
	""" Creates the outputfile, which is a fasta file with all genes from a sp.

	:param sp: current specie.
	:type sp: string
	:param dicoGene: all genes and their sequences, coordinates, etc.
	:type dicoGene: dictionary
	"""
	words = sp.split("_")
	letters = [word[0] for word in words]
	ini = "".join(letters)
	ini = ini.upper()
	output = open("/home/anais/Documents/Data/Genomes/" + sp + "/" + ini + \
		"_gene_unspliced.txt","w")
	for gene in dicoGene :
		chromosome = dicoGene[gene]["Chromosome"]
		start = dicoGene[gene]["geneStart"]
		end = dicoGene[gene]["geneEnd"]
		Sequence = dicoGene[gene]["Sequence"]
		strand = str(dicoGene[gene]["Strand"])
		header = '>'+ gene +' ENSG00000000457 chromosome:'+ \
			dicoGene[gene]['Assembly'] +':'+ chromosome \
			+':'+ str(start) +':'+ str(end) +':'+ strand
		output.write(header + "\n")
		nbLine = math.ceil( float( end - start ) / 60 )
		cpt1 = 0
		cpt2 = 60
		if strand == "-1" :
			Sequence = reverseSequence(Sequence)
		for i in range(0,int(nbLine)) :
			output.write( Sequence[cpt1:cpt2] + "\n" )
			# to have a new line after 60 characters
			cpt1 += 60
			cpt2 += 60
	output.close()

def getJunctionSequences(dicoChromosome, filename):
	dicoJunction = {}
	dicoTr, dicoGene = pgtf.importGTF(filename)
	del dicoGene
	dicoTr = pgtf.getIntron(dicoTr)
	for tr in dicoTr:
		if 'Intron' in dicoTr[tr]:
			for intron in dicoTr[tr]['Intron']:
				if intron != 'Chromosome' and intron != 'Strand':
					if dicoTr[tr]['Chromosome'] in dicoChromosome:
						header = '>'+ tr +' ENSG00000000457 chromosome:'+ \
						dicoTr[tr]['Assembly'] \
						+':'+ dicoTr[tr]['Chromosome'] \
						+':'+ str(dicoTr[tr]['Intron'][intron]['Start']) \
						+':'+ str(dicoTr[tr]['Intron'][intron]['End']) \
						+':'+ dicoTr[tr]['Strand']
						seq = dicoChromosome[dicoTr[tr]['Chromosome']]\
								[dicoTr[tr]['End']-100:dicoTr[tr]['End']]
						dicoJunction[header] = seq
	return dicoJunction

def writeFastaJunction(dicoJunction, sp):
	words = sp.split("_")
	letters = [word[0] for word in words]
	ini = "".join(letters)
	ini = ini.upper()
	output = open("/home/anais/Documents/Data/Genomes/" + sp + "/" + ini + \
		"_transcript_unspliced.txt","w")
	for junction in dicoJunction:
		start = junction.split(':')[3]
		end = junction.split(':')[4]
		strand = junction.split(':')[5]
		sequence = dicoJunction[junction]
		output.write(junction + "\n")
		nbLine = math.ceil( float( len(sequence) ) / 60 )
		cpt1 = 0
		cpt2 = 60
		if strand == "-1" :
			Sequence = reverseSequence(sequence)
		for i in range(0,int(nbLine)) :
			output.write( sequence[cpt1:cpt2] + "\n" )
			# to have a new line after 60 characters
			cpt1 += 60
			cpt2 += 60
	output.close()

def main(sp):
	print "Fasta for " + sp
	dicoChromosome = importFastaChromosome(sp)
	dicoGene = getGeneSequence(sp, dicoChromosome)
	createFastaGene(sp, dicoGene)
	dicoJunction = getJunctionSequences(dicoChromosome, '/home/anais/Documents/Data/Genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae.gtf')
	writeFastaJunction(dicoJunction, sp)
	print "\tDone"

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp = arg.specie # specie to analyse
	main(sp)
