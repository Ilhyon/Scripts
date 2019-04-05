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
	words = sp.split("_")
	letters = [word[0] for word in words]
	ini = "".join(letters)
	ini = ini.upper()
	output = open("/home/anais/Documents/Data/Genomes/" + sp + "/" + ini + \
		"_transcript_unspliced.txt","w") # file opening for reading
	for transcript in dicoFeature :
		nbExon = len(dicoFeature[transcript]["Exon"])
		gene = dicoFeature[transcript]["Gene"]
		for exon in dicoFeature[transcript]["Exon"] :
			chromosome = dicoFeature[transcript]["Exon"][exon]["Chromosome"]
			biotype = dicoFeature[transcript]["Exon"][exon]["Biotype"]
			start = dicoFeature[transcript]["Exon"][exon]["Start"]
			end = dicoFeature[transcript]["Exon"][exon]["End"]
			strand = str(dicoFeature[transcript]["Exon"][exon]["Strand"])
			rank = dicoFeature[transcript]["Exon"][exon]["Rank"]
			start5UTR, end5UTR, start3UTR, end3UTR = \
				retrieveUTRFromDico(dicoFeature[transcript])
			line = gene+"\t"+transcript+"\t"+chromosome+"\t"+biotype+"\t"
			if int(rank) == 1 and int(rank) == nbExon :
				# if the transcript contain only one exon we put all the UTR
				output.write(line
							+start5UTR+"\t"
							+end5UTR+"\t"
							+start3UTR+"\t"
							+end3UTR+"\t"
							+start+"\t"
							+end+"\t"
							+rank+"\t"
							+strand+"\n")
			elif int(rank) == 1 :
				# if it's not the only one exon but it's the first one,
				# we add only the 5UTR
				output.write(line
							+start5UTR+"\t"
							+end5UTR+"\t\t\t"
							+start+"\t"
							+end+"\t"
							+rank+"\t"
							+strand+"\n")
			elif int(rank) == nbExon :
				# if it's not the only one exon but it's the last one,
				# we add only the 3UTR
				output.write(line
							+"\t\t"
							+start3UTR+"\t"
							+end3UTR+"\t"
							+start+"\t"
							+end+"\t"
							+rank+"\t"
							+strand+"\n")
			else :
				# it's only a "midle" exon so we add no UTR
				output.write(line
							+"\t\t\t\t"
							+start+"\t"
							+end+"\t"
							+rank+"\t"
							+strand+"\n")
	output.close()

def writeIDGene(sp, dicoFeature) :
	"""
		Create a file with all gene ID. I don't use this file anymore
		but I let the function just in case.
	"""
	output = open("/home/anais/Documents/Data/Genomes/"+sp+"/"+ini+"_GeneID.txt","w") # file opening for reading
	for gene in dicoFeature :
		output.write(gene+"\n")
	output.close()

def changeStrandFormat(strand) :
	if strand == "+":
		strand = 1
	elif strand == "-":
		strand = -1
	return strand

def retrieveBiotypeFronAttributes(attributes, feature) :
	biotype = ""
	for attribute in attributes:
		if re.search(feature+"_biotype", attribute):
			biotype = attribute.split('"')[1]
	return biotype

def retrieveIdTrFronAttributes(attributes) :
	idTr = ""
	for attribute in attributes:
		if re.search("transcript_id", attribute):
			idTr = attribute.split('"')[1]
	return idTr

def retrieveInfoExonFronAttributes(attributes) :
	rank = ""
	idExon = 0
	for attribute in attributes:
		if re.search("exon_number", attribute):
			rank = attribute.split('"')[1]
		elif re.search("exon_id", attribute):
			idExon = attribute.split('"')[1]
	return rank, idExon

def retrieveUTRFromDico(dico) :
	"""
		From the dictionary of feature, we retrieve the coordinates of
		UTR if they exist
	"""
	start5UTR = ""
	end5UTR = ""
	start3UTR = ""
	end3UTR = ""
	if "5UTR" in dico :
		start5UTR, end5UTR = retrieveCoordUTRFromDico(dico["5UTR"])
	if "3UTR" in dico :
		start3UTR, end3UTR = retrieveCoordUTRFromDico(dico["3UTR"])
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
					startFeature = words[3]
					endFeature = words[4]
					strand = words[6]
					if strand == "+":
						strand = "1"
					elif strand == "-":
						strand = "-1"
					dicoGene[idGene] = {'Chromosome' : chrm,
										'geneStart' : int(startFeature),
										'geneEnd' : int(endFeature),
										'Strand' : strand}
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
					words=l.split('\t')
					attributes = words[8].split(';')
					idGene = attributes[0].split('"')[1]
					feature = words[2]
					chrm = words[0]
					startFeature = words[3]
					endFeature = words[4]
					strand = words[6]
					strand = changeStrandFormat(strand)
					biotype = retrieveBiotypeFronAttributes(attributes, feature)
					if feature == "exon" :
						idTr = retrieveIdTrFronAttributes(attributes)
						rank, idExon = retrieveInfoExonFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, "Exon"))
						dicoTr[idTr]["Exon"].update({idExon : {"Chromosome" : chrm,
																"Start" : startFeature,
																"End" :endFeature,
																"Biotype" : biotype,
																"Strand" : strand,
																"Rank" : rank}})
					elif feature == "five_prime_utr" :
						idTr = retrieveIdTrFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, "5UTR"))
						dicoTr[idTr]["5UTR"].update({"Chromosome" : chrm,
													"Start" : startFeature,
													"End" :endFeature,
													"Biotype" : biotype,
													"Strand" : strand})
					elif feature == "three_prime_utr" :
						idTr = retrieveIdTrFronAttributes(attributes)
						dicoTr.update(createKeyTranscript(dicoTr, idGene, idTr, "3UTR"))
						dicoTr[idTr]["3UTR"].update({"Chromosome" : chrm,
													"Start" : startFeature,
													"End" :endFeature,
													"Biotype" : biotype,
													"Strand" : strand})
					elif feature == "gene":
						dicoGene.update({idGene : {"Chromosome" : chrm,
													"Start" : startFeature,
													"End" :endFeature,
													"Biotype" : biotype,
													"Strand" : strand}})
					elif feature == "transcript":
						idTr = retrieveIdTrFronAttributes(attributes)
						if idTr not in dicoTr:
							dicoTr.update({idTr : {}})
						dicoTr[idTr].update({"Gene" : idGene,
											"Chromosome" : chrm,
											"Start" : startFeature,
											"End" : endFeature,
											"Biotype" : biotype,
											"Strand" : strand})
	else:
		print "This file don't exist : " + filename
	return dicoTr, dicoGene

# if __name__ == '__main__':
# 	parser = build_arg_parser()
# 	arg = parser.parse_args()
# 	sp = arg.specie	# specie to parse
# 	dicoTr, dicoGene = importGTF(sp)
# 	pprint(dicoTr)
