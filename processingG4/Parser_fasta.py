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
import random
import argparse
from Bio import SeqIO
from pprint import pprint
import Parser_gtf as pgtf
import recurrentFunction as rF

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_Fasta')
	parser.add_argument ('-sp', '--specie',
		default = 'yersinia_pestis_biovar_microtus_str_91001')
	return parser

def importGeneFasta(filename, path):
	fastaOrigin = SeqIO.parse(open(filename),'fasta')
	fastaRandom = {}
	for fasta in fastaOrigin:
		name, sequence = fasta.id, str(fasta.seq)
		sequence = list(sequence)
		random.shuffle(sequence)
		fastaRandom[name] = ''.join(sequence)
	output = open(path + 'Random.fas', "w")
	for id in fastaRandom:
		output.write(">" + id + "\n")
		nbLine = math.ceil( float( len(fastaRandom[id]) ) / 60 )
		cpt1 = 0
		cpt2 = 60
		for i in range(0,int(nbLine)) :
			output.write( fastaRandom[id][cpt1:cpt2] + "\n" )
			# to have a new line after 60 characters
			cpt1 += 60
			cpt2 += 60
	output.close()

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
	are those that are present in pan-compara. We also compute the GC content
	of all chromosomes (genome) and print it.

	:param sp: current specie.
	:type sp: string

	:returns: dicoChromosome, all chromosome and their sequences.
	:rtype: dictionary
	"""
	directory = "/home/anais/Documents/Data/Genomes/" + sp + "/Fasta/"
	dicoChromosome = {}
	count = { 'G' : 0, 'C' : 0, 'A' : 0, 'T' : 0}
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
				for nucleotide in count:
					count[nucleotide] += sequence.count(nucleotide)
	print("\tDictionnary fasta done")
	GCcontent = float(count['G'] + count['C']) / \
		(count['G'] + count['C'] + count['A'] + count['T']) * 100
	print('\tGC content in genome is: ' + str(GCcontent) + '.')
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
		dicoGene = pgtf.importGTFGene(filename)
		for geneId in dicoGene:
			chrm = dicoGene[geneId]['Chromosome']
			if chrm in dicoChromosome :
				if (dicoGene[geneId]['Start'] > 0 and
					dicoGene[geneId]['End'] > 0):
					"""
					-1 for coords because in python : start at 0 and not 1
					+1 for the sequence :
					1 2 3 4 5 6 7 8 9
					A T C G G G T A G
					exemple : 6-3 = 3 but there is 4 nucleotides
					so 6-3 +1 = 4
					"""
					geneSequence = \
						dicoChromosome[chrm]\
						[ dicoGene[geneId]['Start']-1:\
						(dicoGene[geneId]['End']-1+1) ]
				elif dicoGene[geneId]['Start'] < 0 :
					# genes overlaping the origin of
					# replication in bacteria
					negSeq = dicoChromosome[chrm] \
						[- -int(dicoGene[geneId]['Start']):]
						# sequence before the origin of replication
					posSeq = dicoChromosome[chrm]\
						[0:int(dicoGene[geneId]['End'])]
						# sequence after the origin of replication
					geneSequence = negSeq + posSeq
				dicoGene[geneId].update({"Sequence" : geneSequence})
	else:
		print('gtf do not exist')
	return dicoGene

def createFasta(sp, dicoFeature, featureType):
	"""Writes a fasta file depending on the feature.

	Fasta file are created for genes or junction. Header are formated like those
	of Ensembl. The identifier is the gene/transcript id. G4RNA Screener uses a
	regex which detects only Human gene/transcript id as stable id. For this
	reasin, the stable id is always ENSG00000000457. If the feature is gene, we
	compute the GC content of all genes and print it.

	:param sp: current specie.
	:type sp: string
	:param dicoFeature: all genes and their sequences, coordinates, etc.
	:type dicoFeature: dictionary
	:param featureType: gene/junction, used to create the correct file dependign on
		the feature
	:type featureType: string
	"""
	ini = rF.setUpperLetter(sp)
	if featureType == 'gene':
		output = open("/home/anais/Documents/Data/Genomes/" + sp + "/" + ini + \
			"_gene_unspliced.txt","w")
		count = { 'G' : 0, 'C' : 0, 'A' : 0, 'T' : 0}
	else:
		output = open("/home/anais/Documents/Data/Genomes/" + sp + "/" + ini + \
			"_transcript_unspliced.txt","w")
	for feature in dicoFeature:
		if 'Sequence' in dicoFeature[feature]:
			if featureType == 'gene':
				chromosome = dicoFeature[feature]["Chromosome"]
				start = dicoFeature[feature]["Start"]
				end = dicoFeature[feature]["End"]
				sequence = dicoFeature[feature]["Sequence"]
				strand = str(dicoFeature[feature]["Strand"])
				header = '>'+ feature +' ENSG00000000457 chromosome:'+ \
					dicoFeature[feature]['Assembly'] +':'+ chromosome \
					+':'+ str(start) +':'+ str(end) +':'+ strand
				for nucleotide in count:
					count[nucleotide] += sequence.count(nucleotide)
			else:
				strand = feature.split(':')[5]
				sequence = dicoFeature[feature]["Sequence"]
				header = feature
			output.write(header + "\n")
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
	try:
	  count
	except:
		pass
	else:
		GCcontent = float(count['G'] + count['C']) / \
			(count['G'] + count['C'] + count['A'] + count['T']) * 100
		print('\tGC content in genes is: ' + str(GCcontent) + '.')

def getJunctionSequences(dicoChromosome, filename):
	"""Gets junction sequence.

	We retrieves the 100 nucleotids upstream an intron and the 100 nucleotides
	downstream of the same intron. This sequence will be a junction sequence.
	Its coordinates will be the start and the end of the intron which is
	spliced.

	:param dicoChromosome: contains all chromosomes sequences.
	:type dicoChromosome: dictionary
	:param filename: name of the gtf file. it will be imported
	:type filename: string

	:returns: dicoJunction, contains all junction and its sequence.
	:rtype: dictionary:
	"""
	dicoJunction = {}
	dicoTr = pgtf.importGTF(filename)
	dicoTr = pgtf.getIntron(dicoTr)
	for tr in dicoTr:
		if 'Intron' in dicoTr[tr]:
			for intron in dicoTr[tr]['Intron']:
				if intron != 'Chromosome' and intron != 'Strand':
					if dicoTr[tr]['Chromosome'] in dicoChromosome:
						tmp = dicoTr[tr]['Intron'][intron]
						chrm = dicoTr[tr]['Chromosome']
						header = '>chromosome:' +  dicoTr[tr]['Assembly'] \
						+':'+ chrm \
						+':'+ str(dicoTr[tr]['Intron'][intron]['Start']) \
						+':'+ str(dicoTr[tr]['Intron'][intron]['End']) \
						+':'+ dicoTr[tr]['Strand']
						seq = dicoChromosome[chrm]\
								[tmp['Start']-100:tmp['Start']]
						seq += dicoChromosome[chrm]\
								[tmp['End']:tmp['End']+100]
						if dicoTr[tr]['Intron']['Strand'] == '-1':
							seq = reverseSequence(seq)
						dicoJunction[header] = {'Sequence' : seq}
	return dicoJunction

def main(sp):
	print("Fasta for " + sp)
	filename = '/home/anais/Documents/Data/Genomes/' + sp + '/' + sp + '.gtf'
	dicoChromosome = importFastaChromosome(sp)
	dicoGene = getGeneSequence(sp, dicoChromosome)
	createFasta(sp, dicoGene, 'gene')
	dicoJunction = getJunctionSequences(dicoChromosome, filename)
	createFasta(sp, dicoJunction, 'junction')
	print("\tDone")

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp = arg.specie # specie to analyse
	ini = rF.setUpperLetter(sp)
	path = '/home/anais/Documents/Data/Genomes/' + sp + '/'
	filename = path + ini + '_gene_unspliced.txt'
	importGeneFasta(filename, path)
