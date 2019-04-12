#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import argparse
import os
from pprint import pprint

# \file Parser_gtf.py
# \author Anaïs Vannutelli 18129080
# \brief This script parse a gtf file to only retrieve informations
# will use later

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_gtf')
	parser.add_argument ('-sp', '--specie', default = 'MM')
	return parser

def writeTranscriptFile(sp, dicoFeature) :
	words = sp.split('_')
	letters = [word[0] for word in words]
	ini = ''.join(letters)
	ini = ini.upper()
	output = open('/home/anais/Documents/Data/Genomes/' + sp + '/' + ini + \
		'_transcript_unspliced.txt','w') # file opening for reading
	for transcript in dicoFeature :
		nbExon = len(dicoFeature[transcript]['Exon'])
		gene = dicoFeature[transcript]['Gene']
		for exon in dicoFeature[transcript]['Exon'] :
			chromosome = dicoFeature[transcript]['Exon'][exon]['Chromosome']
			biotype = dicoFeature[transcript]['Exon'][exon]['Biotype']
			start = dicoFeature[transcript]['Exon'][exon]['Start']
			end = dicoFeature[transcript]['Exon'][exon]['End']
			strand = str(dicoFeature[transcript]['Exon'][exon]['Strand'])
			rank = dicoFeature[transcript]['Exon'][exon]['Rank']
			start5UTR, end5UTR, start3UTR, end3UTR = \
				retrieveUTRFromDico(dicoFeature[transcript])
			line = gene+'\t'+transcript+'\t'+chromosome+'\t'+biotype+'\t'
			if int(rank) == 1 and int(rank) == nbExon :
				# if the transcript contain only one exon we put all the UTR
				output.write(line
							+start5UTR+'\t'
							+end5UTR+'\t'
							+start3UTR+'\t'
							+end3UTR+'\t'
							+start+'\t'
							+end+'\t'
							+rank+'\t'
							+strand+'\n')
			elif int(rank) == 1 :
				# if it's not the only one exon but it's the first one,
				# we add only the 5UTR
				output.write(line
							+start5UTR+'\t'
							+end5UTR+'\t\t\t'
							+start+'\t'
							+end+'\t'
							+rank+'\t'
							+strand+'\n')
			elif int(rank) == nbExon :
				# if it's not the only one exon but it's the last one,
				# we add only the 3UTR
				output.write(line
							+'\t\t'
							+start3UTR+'\t'
							+end3UTR+'\t'
							+start+'\t'
							+end+'\t'
							+rank+'\t'
							+strand+'\n')
			else :
				# it's only a 'midle' exon so we add no UTR
				output.write(line
							+'\t\t\t\t'
							+start+'\t'
							+end+'\t'
							+rank+'\t'
							+strand+'\n')
	output.close()

def writeIDGene(sp, dicoFeature) :
	"""
		Create a file with all gene ID. I don't use this file anymore
		but I let the function just in case.
	"""
	output = open('/home/anais/Documents/Data/Genomes/'+sp+'/'+ini+'_GeneID.txt','w') # file opening for reading
	for gene in dicoFeature :
		output.write(gene+'\n')
	output.close()

def changeStrandFormat(strand) :
	if strand == '+':
		strand = '1'
	elif strand == '-':
		strand = '-1'
	return strand

def retrieveBiotypeFronAttributes(attributes, feature) :
	biotype = ''
	for attribute in attributes:
		if re.search(feature+"_biotype", attribute):
			biotype = attribute.split('"')[1]
	return biotype

def retrieveIdTrFronAttributes(attributes) :
	idTr = ''
	for attribute in attributes:
		if re.search('transcript_id', attribute):
			idTr = attribute.split('"')[1]
	return idTr

def retrieveInfoExonFronAttributes(attributes) :
	rank = 0
	idExon = 0
	for attribute in attributes:
		if re.search('exon_number', attribute):
			rank = int(attribute.split('"')[1])
		elif re.search('exon_id', attribute):
			idExon = attribute.split('"')[1]
	return rank, idExon

def retrieveUTRFromDico(dico) :
	"""
		From the dictionary of feature, we retrieve the coordinates of
		UTR if they exist
	"""
	start5UTR = ''
	end5UTR = ''
	start3UTR = ''
	end3UTR = ''
	if '5UTR' in dico :
		start5UTR, end5UTR = retrieveCoordUTRFromDico(dico['5UTR'])
	if '3UTR' in dico :
		start3UTR, end3UTR = retrieveCoordUTRFromDico(dico['3UTR'])
	return start5UTR, end5UTR, start3UTR, end3UTR

def retrieveCoordUTRFromDico(dico) :
	start = dico["Start"]
	end = dico["End"]
	return start, end

def createKeyTranscript(dico, idGene, idTr, feature) :
	if idTr not in dico :
		dico[idTr] = {"Gene" : idGene,
						feature : {}}
	if feature not in dico[idTr] :
		dico[idTr].update({feature : {}})
	return dico

def importGTFSequence(filename):
	dicoGene = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: # browse all lines
			if not l.startswith('#') and l:
				words = l.split('\t')
				idGene = words[8].split(';')[0].split('"')[1]
				feature = words[2]
				if feature == "gene" :
					chrm = words[0]
					startFeature = int(words[3])
					endFeature = int(words[4])
					strand = words[6]
					if strand == '+':
						strand = '1'
					elif strand == '-':
						strand = '-1'
					dicoGene[idGene] = {'Chromosome' : chrm,
										'geneStart' : startFeature,
										'geneEnd' : endFeature,
										'Strand' : strand,
										'Assembly' : assembly}
			elif re.search('genome-version', l):
				assembly = l.split(' ')[1]
	return dicoGene

def importGTF(filename):
	exists = os.path.isfile(filename)
	if exists :
		dicoTr = {}
		dicoGene = {}
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					words = l.split('\t')
					attributes = words[8].split(';')
					idGene = attributes[0].split('"')[1]
					feature = words[2]
					chrm = words[0]
					startFeature = int(words[3])
					endFeature = int(words[4])
					strand = words[6]
					strand = changeStrandFormat(strand)
					biotype = retrieveBiotypeFronAttributes(attributes, feature)
					if feature == 'exon' :
						idTr = retrieveIdTrFronAttributes(attributes)
						rank, idExon = retrieveInfoExonFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, 'Exon'))
						dicoTr[idTr]['Exon'].update({idExon : {'Chromosome' : chrm,
																'Start' : startFeature,
																'End' : endFeature,
																'Biotype' : biotype,
																'Strand' : strand,
																'Rank' : rank}})
					elif feature == 'five_prime_utr' :
						idTr = retrieveIdTrFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, '5UTR'))
						dicoTr[idTr]['5UTR'].update({'Chromosome' : chrm,
													'Start' : startFeature,
													'End' :endFeature,
													'Biotype' : biotype,
													'Strand' : strand})
					elif feature == 'three_prime_utr' :
						idTr = retrieveIdTrFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, '3UTR'))
						dicoTr[idTr]['3UTR'].update({'Chromosome' : chrm,
													'Start' : startFeature,
													'End' :endFeature,
													'Biotype' : biotype,
													'Strand' : strand})
					elif feature == 'gene':
						dicoGene.update({idGene : {'Chromosome' : chrm,
													'Start' : startFeature,
													'End' :endFeature,
													'Biotype' : biotype,
													'Strand' : strand,
													'Assembly' : assembly}})
					elif feature == 'transcript':
						idTr = retrieveIdTrFronAttributes(attributes)
						if idTr not in dicoTr:
							dicoTr.update({idTr : {}})
						dicoTr[idTr].update({'Gene' : idGene,
											'Chromosome' : chrm,
											'Start' : startFeature,
											'End' : endFeature,
											'Biotype' : biotype,
											'Strand' : strand,
											'Assembly' : assembly})
				elif re.search('genome-version', l):
					assembly = l.split(' ')[1]
	else:
		print "This file don't exist : " + filename
	return dicoTr, dicoGene

def computesCoordRank1(strand, start, end):
	"""Computes the coordinates of the first intron of a transcript.

	Depending of the strand, the start of the first intron is computed. It
	corresponds to the nucleotids next to the first exon end.

	:param strand: strand of a transcript.
	:type strand: string
	:param start: chromosomal start of the exon.
	:type start: integer
	:param end: chromosomal end of the exon.
	:type end: integer

	:returns: intronStart, start of the first exon
	:rtype: integer
	"""
	if strand == '1':
		intronStart = end + 1
	else:
		intronStart = start - 1
	return intronStart

def computesCoordLastRank(strand, start, end):
	"""Computes the coordinates of the last intron of a transcript.

	Depending of the strand, the end of the last intron is computed. It
	corresponds to the nucleotids before the last exon start.

	:param strand: strand of a transcript.
	:type strand: string
	:param start: chromosomal start of the exon.
	:type start: integer
	:param end: chromosomal end of the exon.
	:type end: integer

	:returns: intronEnd, end of the last exon
	:rtype: integer
	"""
	if strand == '1':
		intronEnd = start - 1
	else:
		intronEnd = end + 1
	return intronEnd

def computesCoordOther(strand, start, end):
	"""Computes the coordinates of midle intron of a transcript.

	With the coordinates of a "midle" exon we can compute the end of the
	previous exon and then the start of the next intron.

	:param strand: strand of a transcript.
	:type strand: string
	:param start: chromosomal start of the exon.
	:type start: integer
	:param end: chromosomal end of the exon.
	:type end: integer

	:returns: intronStart and intronEnd, end of the previous intron and start
		of the next intron.
	:rtype: integer
	"""
	if strand == '1':
		intronEnd = start - 1
		intronStart = end + 1
	else:
		intronEnd = end + 1
		intronStart = start - 1
	return intronStart, intronEnd

def setCoordReverseStrand(dico):
	"""Change the coords depending on the strand.

	I choosed to get chromosomal coords on my file. It will always be the
	smaller coords before the superior one. So for a reverse feature
	the start and the end are inversed.

	:param dicoTr: contains all transcript and its feature.
	:type dicoTr: dictionary

	:returns: dicoTr, updated with good coordinates.
	:rtype: dictionary
	"""
	for intron in dico:
		if intron != 'Chromosome' and intron != 'Strand' :
			if (dico['Strand'] == '-1' and
				dico[intron]['End'] < dico[intron]['Start']):
				tmp = dico[intron]['Start']
				dico[intron]['Start'] = dico[intron]['End']
				dico[intron]['End'] = tmp
	return dico

def getIntron(dicoTr):
	for tr in dicoTr:
		if tr != 'Assembly':
			if len(dicoTr[tr]['Exon']) > 2:
				# if there is more then one exon, we compute
				# intron coordinates.
				if 'Intron' not in dicoTr[tr]:
					dicoTr[tr]['Intron'] = {}
				for exon in dicoTr[tr]['Exon']:
					tmp = dicoTr[tr]['Exon'][exon]
					if tmp['Rank'] == 1:
						if tmp['Rank'] not in dicoTr[tr]['Intron']:
							dicoTr[tr]['Intron'][tmp['Rank']] = {}
						dicoTr[tr]['Intron'][1].update( \
							{'Start' : computesCoordRank1(tmp['Strand'],
							tmp['Start'], tmp['End'])})
					elif tmp['Rank'] == len(dicoTr[tr]['Exon']):
						if tmp['Rank']-1 not in dicoTr[tr]['Intron']:
							dicoTr[tr]['Intron'][tmp['Rank']-1] = {}
						dicoTr[tr]['Intron'][tmp['Rank']-1].update( \
							{'End' : computesCoordLastRank(tmp['Strand'],
							tmp['Start'], tmp['End'])})
					else:
						if tmp['Rank']-1 not in dicoTr[tr]['Intron']:
							dicoTr[tr]['Intron'][ tmp['Rank']-1 ] = {}
						if tmp['Rank'] not in dicoTr[tr]['Intron']:
							dicoTr[tr]['Intron'][tmp['Rank']] = {}
						intronStart, intronEnd = \
							computesCoordOther(tmp['Strand'], \
							tmp['Start'], tmp['End'])
						dicoTr[tr]['Intron'][ tmp['Rank'] ].update(\
							{'Start' : intronStart})
						dicoTr[tr]['Intron'][ tmp['Rank']-1 ].update(\
							{'End' : intronEnd})
					dicoTr[tr]['Intron']['Chromosome'] = tmp['Chromosome']
					dicoTr[tr]['Intron']['Strand'] = tmp['Strand']
				dicoTr[tr]['Intron'].update(\
					setCoordReverseStrand(dicoTr[tr]['Intron']))
			elif len(dicoTr[tr]['Exon']) > 1:
				if 'Intron' not in dicoTr[tr]:
					dicoTr[tr]['Intron'] = {}
				for exon in dicoTr[tr]['Exon']:
					tmp = dicoTr[tr]['Exon'][exon]
					if 1 not in dicoTr[tr]['Intron']:
						dicoTr[tr]['Intron'][1] = {}
					if tmp['Rank'] == 1:
						dicoTr[tr]['Intron'][1].update( \
							{'Start' : computesCoordRank1(tmp['Strand'],
							tmp['Start'], tmp['End'])})
					else:
						dicoTr[tr]['Intron'][1].update( \
							{'End' : computesCoordLastRank(tmp['Strand'],
							tmp['Start'], tmp['End'])})
					dicoTr[tr]['Intron']['Chromosome'] = \
						tmp['Chromosome']
					dicoTr[tr]['Intron']['Strand'] = tmp['Strand']
				dicoTr[tr]['Intron'].update(\
					setCoordReverseStrand(dicoTr[tr]['Intron']))
	return dicoTr

def countIntron(dicoTr):
	lenIntron = []
	for tr in dicoTr:
		if 'Intron' in dicoTr[tr]:
			for intron in dicoTr[tr]['Intron']:
				if intron != 'Chromosome' and intron != 'Strand':
					lenIntron.append(dicoTr[tr]['Intron'][intron]['End'] - \
						dicoTr[tr]['Intron'][intron]['Start'])
	print lenIntron

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp = arg.specie	# specie to parse
	filename = "/home/anais/Documents/Data/Genomes/" + sp + \
		"/" + sp + ".gtf"
	dicoTr, dicoGene = importGTF(filename)
	dicoTr = getIntron(dicoTr)
	countIntron(dicoTr)
