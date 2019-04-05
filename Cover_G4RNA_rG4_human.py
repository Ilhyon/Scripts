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

def getStrand(strand):
	if strand == "+":
		strand = "1"
	else :
		strand = "-1"
	return strand

def readLinerG4seq(line):
	words = line.split('\t')
	dicoL = {"Sequence" : words[9],
			"idTr" : words[10],
			"Chromosome" : words[0].split("r")[1],
			"Start" : words[1],
			"End" : words[2],
			"Strand" : words[5],
			"Location" : words[12]}
	dicoL["Strand"] = getStrand(dicoL["Strand"])
	return dicoL
	
def readLineG4RNA(line):
	words = line.split('\t')
	if words[2].split("r")[0] :
		chrm = words[2].split("r")[1]
	else :
		chrm = ""
	dicoL = {"Chromosome" : chrm,
			"geneID" : words[0],
			"g4Type" : words[1],
			"Length" : words[5],
			"Sequence" : words[7],
			"Folding" : words[8],
			"Start" : words[3],
			"End" : words[4],
			"Strand" : words[6]}
	dicoL["Strand"] = getStrand(dicoL["Strand"])
	return dicoL

def readQuery(line):
	words = line.split('\t')
	header = words[0]
	l = {"Sequence" : words[3],
		"idTr" : header.split("|")[0],
		"pG4Id" : header.split("|")[1],
		"Strand" : header.split("|")[2],
		"Biotype" : words[6],
		"Location" : words[5]}
	l.update({"Chromosome" : l["pG4Id"].split(":")[0],
			"Start" : l["pG4Id"].split(":")[1].split("-")[0],
			"End" : l["pG4Id"].split(":")[1].split("-")[1]})
	return l

def removeDoubleG4(dicoG4, listeG4):
	for G4_to_Delete in listeG4 :
		chrm = G4_to_Delete[0]
		coordG4 = G4_to_Delete[1]
		if coordG4 in dicoG4[chrm] :
			del dicoG4[chrm][coordG4]
	return dicoG4

def countG4Target(dicoG4, target):
	cpt = 0
	for chrm in dicoG4 :
		cpt += len(dicoG4[chrm])
	print "Number of G4 in "+ target +" : "+ str(cpt)

def importG4RNA():
	"""
		Read a file that contain all G4 from G4RNA screener (file from 31/01/2019)
		The db maybe update later so the results could change.
		Return a dictionary of G4 folded and WT : 
		{chr : {start-end : {"Start","End", "Sequence"}}
	"""
	directory = "/home/anais/Documents/Data/G4RNA/G4rna_G4.txt"
	dicoG4RNA = {}
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		list_G4_to_Delete = []
		for l in lines:
			if not l.startswith('Gene') and l:
				l = readLineG4RNA(l)
				if (l["Folding"] == "1" and
					l["geneID"] != "Artificial"
					and re.search("WT", l["g4Type"])):
					if l["Chromosome"] not in dicoG4RNA :
						dicoG4RNA[l["Chromosome"]] = {str(l["Start"])+"-"+str(l["End"]) : {}}
					else :
						for G4 in dicoG4RNA[l["Chromosome"]] :
							start_dico = G4.split("-")[0]
							end_dico = G4.split("-")[1]
							if ((l["Start"] == start_dico or l["End"] == end_dico) and
								l["Length"] > dicoG4RNA[l["Chromosome"]][G4]["Length"]):
								list_G4_to_Delete.append([l["Chromosome"], G4])
					dicoG4RNA[l["Chromosome"]].update({str(l["Start"])+"-"+str(l["End"]) : {"Sequence" : l["Sequence"],
													"Start" : l["Start"],
													"End" : l["End"],
													"Length" : l["Length"],
													"Strand" : l["Strand"]}})
	dicoG4RNA = removeDoubleG4(dicoG4RNA, list_G4_to_Delete)
	return dicoG4RNA

def importrG4():
	directory = "/home/anais/Documents/Publis/rG4seq/cdt_K.csv"
	dicorG4 = {}
	with open(directory) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if not l.startswith('chr\t') and l:
				l = readLinerG4seq(l)
				if l["Chromosome"] not in dicorG4 :
					dicorG4[l["Chromosome"]] = {str(l["Start"])+"-"+str(l["End"]) : {}}
				dicorG4[l["Chromosome"]].update({str(l["Start"])+"-"+str(l["End"]) : {"Sequence" : l["Sequence"],
													"Start" : l["Start"],
													"End" : l["End"],
													"Strand" : l["Strand"]}})
	return dicorG4

def importpG4(directoryQuery):
	dicopG4 = {}
	with open(directoryQuery) as f:
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if not l.startswith('Info') and l:
				l = readQuery(l)
				if l["Chromosome"] not in dicopG4 :
					dicopG4[l["Chromosome"]] = {}
				dicopG4[l["Chromosome"]].update({str(l["Start"])+"-"+str(l["End"]) : l})
	return dicopG4

def coverage(dicoTarget, strTarget, dicoQuery, strQuery):
	"""
	Target = bd or rG4
	Query = what I want to compare (pG4)
	"""
	resultTarget = []
	resultQuery = []
	for chrm in dicoQuery:
		if chrm in dicoTarget :
			for queryG4 in dicoQuery[chrm]:
				tmp = dicoQuery[chrm][queryG4]
				if chrm in dicoTarget:
					for targetG4 in dicoTarget[chrm]:
						if tmp["Strand"] == dicoTarget[chrm][targetG4]["Strand"]:
							# now G4 target and query are on the same
							# strand and same chromosome
							coord = targetG4.split("|")[0]
							startTargetG4 = coord.split("-")[0]
							endTargetG4 = coord.split("-")[1]
							targetG4seq = dicoTarget[chrm][targetG4]["Sequence"]
							if targetG4+"\t"+targetG4seq not in resultTarget :
								if ((startTargetG4 > tmp["Start"] and startTargetG4 < tmp["End"]) or
									(endTargetG4 > tmp["Start"] and endTargetG4 < tmp["End"])):
									resultTarget.append(targetG4+"\t"+targetG4seq)
									resultQuery.append(queryG4+"\t"+tmp["Sequence"])
							if targetG4+"\t"+targetG4seq not in resultTarget :
								if ((tmp["Start"] > startTargetG4 and tmp["Start"] < endTargetG4) or
									(tmp["End"] > startTargetG4 and tmp["End"] < endTargetG4)):
									resultTarget.append(targetG4+"\t"+targetG4seq)
									resultQuery.append(queryG4+"\t"+tmp["Sequence"])
							if targetG4+"\t"+targetG4seq not in resultTarget :
								if (re.search(targetG4seq, tmp["Sequence"]) or
									re.search(tmp["Sequence"], targetG4seq)):
									resultTarget.append(targetG4+"\t"+targetG4seq)
									resultQuery.append(queryG4+"\t"+tmp["Sequence"])
	print "Number of G4 from query ("+strQuery+") corresponding to G4 from "+strTarget+" : "+str(len(resultTarget))
	print "Number of G4 from target ("+strTarget+") corresponding to G4 from "+strQuery+" : "+str(len(list(set(resultQuery))))

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_Fasta')
	parser.add_argument ('-p', '--path', default = '/home/anais/Documents/Data/Human/G4Conserve-master/results/all/ini_HS_allScore_G4InTranscript.txt')
	parser.add_argument ('-t', '--type', default = 'allScore')
	parser.add_argument ('-sp', '--specie', default = 'HS')
	return parser

def main(directorypG4, dataType) :
	dicoG4RNA = importG4RNA()
	countG4Target(dicoG4RNA, "G4RNA")
	dicopG4 = importpG4(directorypG4)
	countG4Target(dicopG4, dataType)
	dicorG4 = importrG4()
	countG4Target(dicorG4, "rG4")
	coverage(dicoG4RNA, "G4RNA", dicopG4, dataType)
	coverage(dicoG4RNA, "G4RNA", dicorG4, "rG4")
	coverage(dicorG4, "rG4", dicopG4, dataType)

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	directoryAll = path
	main(directoryAll, "all score")
	
	
