#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import argparse
import os
from pprint import pprint
import recurrentFunction as rF

# \file Parser_gtf.py
# \author Anaïs Vannutelli 18129080
# \brief This script parse a gtf file to only retrieve informations
# will use later

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_gtf')
	parser.add_argument ('-sp', '--specie', default = 'mus_musculus')
	return parser

def changeStrandFormat(strand):
	"""Changes the format of the strand from +/- to 1/-1.
	"""
	if strand == '+':
		strand = '1'
	elif strand == '-':
		strand = '-1'
	return strand

def retrieveBiotypeFronAttributes(feature, attributes):
	"""Gets the biotype from attributes.

	:param feature: transcript, gene, exon, or utr.
	:type feature: string
	:param attributes: last colons of gtf file, contains a lot of informations
		but is different depending on the feature.
	:type attributes: list

	:returns: biotype of the feature. The biotype of a gene can be different
		compared to the transcript biotype.
	:rtype: string
	"""
	biotype = ''
	for attribute in attributes:
		if re.search(feature+"_biotype", attribute):
			biotype = attribute.split('"')[1]
	return biotype

def retrieveIdTrFronAttributes(attributes):
	"""Gets the id transcript from attributes.

	:param attributes: last colons of gtf file, contains a lot of informations
		but is different depending on the feature.
	:type attributes: list

	:returns: idTr, transctipt identifier.
	:rtype: string
	"""
	idTr = ''
	for attribute in attributes:
		if re.search('transcript_id', attribute):
			idTr = attribute.split('"')[1]
	return idTr

def retrieveInfoExonFronAttributes(attributes):
	"""Gets rank and exon number from attributes.

	:param attributes: last colons of gtf file, contains a lot of informations
		but is different depending on the feature.
	:type attributes: list

	:returns: rank, idExon of an exon.
	:rtype: strings
	"""
	rank = 0
	idExon = 0
	for attribute in attributes:
		if re.search('exon_number', attribute):
			rank = int(attribute.split('"')[1])
		elif re.search('exon_id', attribute):
			idExon = attribute.split('"')[1]
	return rank, idExon

def retrieveUTRFromDico(dico):
	"""Gets UTR coords if they exist.

	:param dico: contains all transcripts and their features.
	:type dico: dictionary

	:returns: start5UTR, end5UTR, start3UTR, end3UTR.
	:rtypes: strings
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

def retrieveCoordUTRFromDico(dico):
	start = dico["Start"]
	end = dico["End"]
	return start, end

def createKeyTranscript(dico, idGene, idTr, feature):
	if idTr not in dico :
		dico[idTr] = {"Gene" : idGene,
						feature : {}}
	if feature not in dico[idTr] :
		dico[idTr].update({feature : {}})
	return dico

def importGTFGene(filename):
	"""Imports genes features from the gtf file.

	:param filename: name of the gtf file.
	:type filename: string

	:returns: dicoGene, contains all informations abotu genes.
	:rtype: dictionary
	"""
	dicoGene = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: # browse all lines
			if not l.startswith('#') and l:
				words = l.split('\t')
				idGene = words[8].split(';')[0].split('"')[1]
				attributes = words[8].split(';')
				feature = words[2]
				if feature == "gene" :
					chrm = words[0]
					startFeature = int(words[3])
					endFeature = int(words[4])
					strand = words[6]
					biotype = retrieveBiotypeFronAttributes(feature,
						attributes)
					if strand == '+':
						strand = '1'
					elif strand == '-':
						strand = '-1'
					dicoGene[idGene] = {'Chromosome' : chrm,
										'Start' : startFeature,
										'End' : endFeature,
										'Strand' : strand,
										'Biotype' : biotype,
										'Assembly' : assembly}
			elif re.search('genome-version', l):
				assembly = l.split(' ')[1]
	return dicoGene

def importGTF(filename):
	"""Imports gtf file.

	This function aims to retrieve as more as possible informations contained
	in the gtf file related to transcripts.

	:param filelame: name of the gtf file.
	:type filename: string

	:returns: dicoTr, contains all transcripts from a psecie and all features
		linked to it : UTR, exon, intron.
	:rtype: dictionary
	"""
	exists = os.path.isfile(filename)
	if exists :
		dicoTr = {}
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
					biotype = retrieveBiotypeFronAttributes(feature, attributes)
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
	return dicoTr

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
	"""Gets all intron of a specie transcriptom.

	We browse all transcripts of a specie and computes all intron coordinates.
	Coordinates will be stored as chromosomal coordinates. This mean that the
	start will always be inferior compare to the end.

	:param dicoTr: contains all transcripts from a psecie and all features
		linked to it : UTR, exon, intron.
	:type dicoTr: dictionary

	:returns: dicoTr, updated with introns.
	:rtype: dictionary
	"""
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

def countIntron(dicoTr, sp):
	"""Counts the length of intron and write them in a file.

	Some transcript have intron of 1 nucleotides. We computes the length of
	introns to see for each species what's the range of the intron length.

	:param dicoTr: contains all transcript of a specie and all features in it
		(exon, utr, intron and coords).
	:type dicoTr: dictionnary
	:param sp: name of the specie.
	:type sp: string
	"""
	lenIntron = []
	ini = rF.setUpperLetter(sp)
	for tr in dicoTr:
		if 'Intron' in dicoTr[tr]:
			for intron in dicoTr[tr]['Intron']:
				if intron != 'Chromosome' and intron != 'Strand':
					lenIntron.append(dicoTr[tr]['Intron'][intron]['End'] - \
						dicoTr[tr]['Intron'][intron]['Start'])
	if lenIntron:
		output = open('/home/anais/Documents/Data/Genomes/' + sp + '/' + ini + \
			'_intron_length.txt','w')
		output.write(sp + '\n')
		output.write('\n'.join(str(x) for x in lenIntron))

def writeIntron(sp, dico):
	"""Writes a file with all intron coords.

	:param sp: name of the specie.
	:type sp: string
	:param dico: contains all transcript and its features (utr, exon, intron).
	:type dico: dictionary
	"""
	ini = rF.setUpperLetter(sp)
	output = open('/home/anais/Documents/Data/Genomes/' + sp + '/' + ini + \
		'_intron.txt','w')
	for tr in dico:
		if 'Intron' in dicoTr[tr]:
			for intron in dicoTr[tr]['Intron']:
				if intron != 'Chromosome' and intron != 'Strand':
					line = tr + '\t' + dico[tr]['Chromosome'] + '\t' + \
						str(dico[tr]['Intron'][intron]['Start']) + '\t' + \
						str(dico[tr]['Intron'][intron]['End']) + '\t' + \
						dico[tr]['Strand'] + '\n'
					output.write(line)
	output.close()

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp = arg.specie	# specie to parse
	print sp
	filename = "/home/anais/Documents/Data/Genomes/" + sp + \
		"/" + sp + ".gtf"
	dicoTr, dicoGene = importGTF(filename)
	dicoTr = getIntron(dicoTr)
	writeIntron(sp, dicoTr)
	countIntron(dicoTr, sp)
