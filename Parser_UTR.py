#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import os
import argparse
from pprint import pprint

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_gtf')
	parser.add_argument ('-s', '--specie', default = 'MM')
	return parser
	
def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp=arg.specie	# specie to analyse
	words = sp.split("_")
	letters = [word[0] for word in words]
	ini = "".join(letters)
	ini = ini.upper()
	print 'Specie : '+sp+", Ini : " + ini
	filename = "/home/anais/Documents/Data/Genomes/"+sp+"/"+ini+"_transcript_unspliced.txt"
	exists = os.path.isfile(filename)
	if exists :	
		dicoInfo = {}
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					words = l.split('\t')
					if len(words) != 12:
						print words
					idTr = words[1]
					idGene = words[0]
					chrm = words[2]
					Bt = words[3]
					if words[4] or words[5]:
						fUTRs = words[4]
						fUTRe = words[5]
					else :
						fUTRs = ""
						fUTRe = ""
					if words[6] or words[7]:
						tUTRs = words[6]
						tUTRe = words[7]
					else :
						tUTRs = ""
						tUTRe = ""
					if words[8] or words[9] or words[10] :
						exonS = words[8]
						exonE = words[9]
						rank = words[10]
					else :
						exonS = ""
						exonE = ""
						rank = ""
					strand = words[11]
					if (tUTRs or tUTRe) and (exonS or exonE or rank) and (fUTRs or fUTRe):
						if idTr not in dicoInfo :
							exonName = "Exon " + rank
							dicoInfo.update({idTr : {"cpt" : 1, exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							dicoInfo[idTr].update({"3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							dicoInfo[idTr].update({"5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
						else :
							exonName = "Exon " + rank
							dicoInfo[idTr].update({exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							dicoInfo[idTr]["cpt"] += 1
							dicoInfo[idTr].update({"3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							dicoInfo[idTr].update({"5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
					elif (fUTRs or fUTRe) and (exonS or exonE or rank) :
						if idTr not in dicoInfo:
							exonName = "Exon " + rank
							dicoInfo.update({idTr : {"cpt" : 1, exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							dicoInfo[idTr].update({"5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
						else :
							exonName = "Exon " + rank
							dicoInfo[idTr].update({exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							dicoInfo[idTr]["cpt"] += 1
							dicoInfo[idTr].update({"5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
					elif (tUTRs or tUTRe) and (exonS or exonE or rank):
						if idTr not in dicoInfo :
							exonName = "Exon " + rank
							dicoInfo.update({idTr : {"cpt" : 1, exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							dicoInfo[idTr].update({"3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
						else :
							exonName = "Exon " + rank
							dicoInfo[idTr].update({exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							dicoInfo[idTr]["cpt"] += 1
							dicoInfo[idTr].update({"3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
					else :
						if idTr not in dicoInfo :
							if exonS or exonE or rank :
								exonName = "Exon " + rank
								dicoInfo.update({idTr : {"cpt" : 1, exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							elif tUTRs or tUTRe :
								dicoInfo.update({idTr : {"cpt" : 0, "3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							elif fUTRs or fUTRe :
								dicoInfo.update({idTr : {"cpt" : 0, "5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							else :
								"Ah"
						else :
							if exonS or exonE or rank :
								exonName = "Exon " + rank
								dicoInfo[idTr].update({exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
								dicoInfo[idTr]["cpt"] += 1
							elif tUTRs or tUTRe :
								dicoInfo[idTr].update({"3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							elif fUTRs or fUTRe :
								dicoInfo[idTr].update({"5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							else :
								"Ah"
		words = sp.split("_")
		letters = [word[0] for word in words]
		ini = "".join(letters)
		ini = ini.upper()
		output = open("/home/anais/Documents/Data/Genomes/"+sp+"/"+ini+"_transcript_unspliced.txt","w") # file opening for reading
		for tr in dicoInfo :
			for feature in dicoInfo[tr]:
				if feature != "cpt" and feature != "3UTR" and feature != "5UTR" :
					if "3UTR" in dicoInfo[tr] and "5UTR" in dicoInfo[tr] :
						if feature == "Exon 1" and feature == "Exon "+str(dicoInfo[tr]["cpt"]):
							output.write(dicoInfo[tr]["Exon 1"]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr]["Exon 1"]["Chromosome"]+"\t"+dicoInfo[tr]["Exon 1"]["Biotype"]+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRstart"])+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRend"])+"\t"+str(dicoInfo[tr]["3UTR"]["3UTRstart"])+"\t"+str(dicoInfo[tr]["3UTR"]["3UTRend"])+"\t"+str(dicoInfo[tr]["Exon 1"]["exonS"])+"\t"+str(dicoInfo[tr]["Exon 1"]["exonE"])+"\t"+str(dicoInfo[tr]["Exon 1"]["Rank"])+"\t"+dicoInfo[tr]["Exon 1"]["Strand"]+"\n")
						elif feature == "Exon 1" :
							output.write(dicoInfo[tr]["Exon 1"]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr]["Exon 1"]["Chromosome"]+"\t"+dicoInfo[tr]["Exon 1"]["Biotype"]+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRstart"])+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRend"])+"\t\t\t"+str(dicoInfo[tr]["Exon 1"]["exonS"])+"\t"+str(dicoInfo[tr]["Exon 1"]["exonE"])+"\t"+str(dicoInfo[tr]["Exon 1"]["Rank"])+"\t"+dicoInfo[tr]["Exon 1"]["Strand"]+"\n")
						elif feature == "Exon "+str(dicoInfo[tr]["cpt"]) :
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t"+str(dicoInfo[tr]["3UTR"]["3UTRstart"])+"\t"+str(dicoInfo[tr]["3UTR"]["3UTRend"])+"\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
						elif feature !=  "Exon 1" and feature != "Exon "+str(dicoInfo[tr]["cpt"]):
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t\t\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
						else :
							print "Ah"
					elif "3UTR" in dicoInfo[tr] :
						if feature == "Exon "+str(dicoInfo[tr]["cpt"]) :
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t"+str(dicoInfo[tr]["3UTR"]["3UTRstart"])+"\t"+str(dicoInfo[tr]["3UTR"]["3UTRend"])+"\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
						if feature != "Exon "+str(dicoInfo[tr]["cpt"]):
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t\t\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
					elif "5UTR" in dicoInfo[tr] :
						if feature == "Exon 1" :
							output.write(dicoInfo[tr]["Exon 1"]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr]["Exon 1"]["Chromosome"]+"\t"+dicoInfo[tr]["Exon 1"]["Biotype"]+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRstart"])+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRend"])+"\t\t\t"+str(dicoInfo[tr]["Exon 1"]["exonS"])+"\t"+str(dicoInfo[tr]["Exon 1"]["exonE"])+"\t"+str(dicoInfo[tr]["Exon 1"]["Rank"])+"\t"+dicoInfo[tr]["Exon 1"]["Strand"]+"\n")
						if feature !=  "Exon 1":
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t\t\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
					else:
						output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+str(dicoInfo[tr][feature]["Chromosome"])+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t\t\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+str(dicoInfo[tr][feature]["Strand"])+"\n")
		output.close()
		print "Done"
	else :
		print "Echec for : "+sp
# ~ pprint(dicoInfo)

main()














