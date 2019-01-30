#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import argparse
import os
from pprint import pprint

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Parser_gtf')
	parser.add_argument ('-s', '--specie', default = 'MM')
	return parser
	
def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp=arg.specie	# specie to analyse
	filename = "/home/anais/Documents/Data/Genomes/"+sp+"/"+sp+".gtf"
	exists = os.path.isfile(filename)
	if exists :	
		print "GTF for "+sp
		dicoFeature = {}
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					words=l.split('\t')
					attribute = words[8].split(';')
					idGene = attribute[0].split('"')[1]
					feature = words[2]
					chrm = words[0]
					startFeature = words[3]
					endFeature = words[4]
					strand = words[6]
					if strand == "+":
						strand = 1
					elif strand == "-":
						strand = -1
					biotype = ""
					for i in range(0,len(attribute)):
						if re.search("transcript_biotype", attribute[i]):
							biotype = attribute[i].split('"')[1]
					if idGene not in dicoFeature :
						if feature == "gene" :
							dicoFeature[idGene] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
						elif feature == "transcript" :
							for i in range(0,len(attribute)):
								if re.search("transcript_id", attribute[i]):
									idTr = attribute[i].split('"')[1]
							dicoFeature[idGene] = {"Transcript" : idTr}
							dicoFeature[idGene]["Transcript"][idTr] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
						elif feature == "exon" :
							rank = ""
							idExon = 0
							idTr = 0
							for i in range(0,len(attribute)):
								if re.search("exon_number", attribute[i]):
									rank = attribute[i].split('"')[1]
								elif re.search("exon_id", attribute[i]):
									idExon = attribute[i].split('"')[1]
								elif re.search("transcript_id", attribute[i]):
									idTr = attribute[i].split('"')[1]
							dicoFeature[idGene] = {"Transcript" : idTr}
							dicoFeature[idGene]["Transcript"][idTr] = {"Exon" : idExon}
							dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
						elif feature == "five_prime_utr" :
							idTr = 0
							for i in range(0,len(attribute)):
								if re.search("transcript_id", attribute[i]):
									idTr = attribute[i].split('"')[1]
							dicoFeature[idGene] = {"Transcript" : idTr}
							dicoFeature[idGene]["Transcript"][idTr]={"5UTR" :{"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}}
						elif feature == "three_prime_utr" :
							idTr = 0
							for i in range(0,len(attribute)):
								if re.search("transcript_id", attribute[i]):
									idTr = attribute[i].split('"')[1]
							dicoFeature[idGene] = {"Transcript" : idTr}
							dicoFeature[idGene]["Transcript"][idTr]={"3UTR" :{"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}}
					else:
						if feature == "gene" :
							print "This gene already exist :"
							pprint(dicoFeature[idGene])
							print "Chromosome : "+chrm+" | Start : "+startFeature+" | End :"+endFeature+" | Biotype : "+biotype+" | Strand : "+strand
						elif feature == "transcript" :
							idTr = 0
							for i in range(0,len(attribute)):
								if re.search("transcript_id", attribute[i]):
									idTr = attribute[i].split('"')[1]
							if "Transcript" in dicoFeature[idGene]:
								if idTr not in dicoFeature[idGene]["Transcript"] :
									dicoFeature[idGene]["Transcript"][idTr] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
								else :
									print "This transcript already exist :"
									print idTr
									pprint(dicoFeature[idGene]["Transcript"][idTr])
									print "Chromosome : "+chrm+" | Start : "+startFeature+" | End :"+endFeature+" | Biotype : "+biotype+" | Strand : "+strand
							else : 
								dicoFeature[idGene]["Transcript"] = {idTr:{}}
								dicoFeature[idGene]["Transcript"][idTr] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
						elif feature == "exon" :
							rank = ""
							idExon = 0
							idTr = 0
							for i in range(0,len(attribute)):
								if re.search("exon_number", attribute[i]):
									rank = attribute[i].split('"')[1]
								elif re.search("exon_id", attribute[i]):
									idExon = attribute[i].split('"')[1]
								elif re.search("transcript_id", attribute[i]):
									idTr = attribute[i].split('"')[1]
							if "Transcript" in dicoFeature[idGene]:
								if idTr not in dicoFeature[idGene]["Transcript"] :
									dicoFeature[idGene]["Transcript"].update({idTr : {}})
									dicoFeature[idGene]["Transcript"][idTr]["Exon"] = {idExon:{}}
									dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
								else :
									if "Exon" in dicoFeature[idGene]["Transcript"][idTr]:
										if idExon not in dicoFeature[idGene]["Transcript"][idTr]["Exon"]:
											dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
										else:
											print "This exon already exist :"
											pprint(dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon])
											print "Chromosome : "+chrm+" | Start : "+startFeature+" | End :"+endFeature+" | Biotype : "+biotype+" | Strand : "+strand
									else :
										dicoFeature[idGene]["Transcript"][idTr]["Exon"] ={idExon : {}}
										dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
							else : 
								dicoFeature[idGene].update({"Transcript" : idTr})
								dicoFeature[idGene]["Transcript"][idTr] = {"Exon" : idExon}
								dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
						elif feature == "five_prime_utr" :
							idTr = 0
							for i in range(0,len(attribute)):
								if re.search("transcript_id", attribute[i]):
									idTr = attribute[i].split('"')[1]
							if "Transcript" in dicoFeature[idGene]:
								if idTr not in dicoFeature[idGene]["Transcript"] :
									dicoFeature[idGene]["Transcript"].update({idTr : {"5UTR" :{"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}}})
								else :
									if "5UTR" in dicoFeature[idGene]["Transcript"][idTr] :
										if (strand == "1" and startFeature == dicoFeature[idGene]["Transcript"][idTr]["Start"]) or (strand == "-1" and endFeature == dicoFeature[idGene]["Transcript"][idTr]["End"]):
											dicoFeature[idGene]["Transcript"][idTr]["5UTR"] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
									else :
										dicoFeature[idGene]["Transcript"][idTr]["5UTR"] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
							else : 
								dicoFeature[idGene].update({"Transcript" : idTr})
								dicoFeature[idGene]["Transcript"][idTr] = {"5UTR" : {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}}
						elif feature == "three_prime_utr" :
							idTr = 0
							for i in range(0,len(attribute)):
								if re.search("transcript_id", attribute[i]):
									idTr = attribute[i].split('"')[1]
							if "Transcript" in dicoFeature[idGene]:
								if idTr not in dicoFeature[idGene]["Transcript"] :
									dicoFeature[idGene]["Transcript"].update({idTr : {"3UTR" :{"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}}})
								else :
									if "3UTR" in dicoFeature[idGene]["Transcript"][idTr]:
										if (strand == "1" and endFeature == dicoFeature[idGene]["Transcript"][idTr]["End"]) or (strand == "-1" and startFeature == dicoFeature[idGene]["Transcript"][idTr]["Start"]):
											dicoFeature[idGene]["Transcript"][idTr]["3UTR"] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
									else :
										dicoFeature[idGene]["Transcript"][idTr]["3UTR"] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
							else : 
								dicoFeature[idGene].update({"Transcript" : idTr})
								dicoFeature[idGene]["Transcript"][idTr] = {"3UTR" : {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}}
		
		#~ pprint(dicoFeature)
		words = sp.split("_")
		letters = [word[0] for word in words]
		ini = "".join(letters)
		ini = ini.upper()
		output = open("/home/anais/Documents/Data/Genomes/"+sp+"/"+ini+"_transcript_unspliced.txt","w") # file opening for reading
		for gene in dicoFeature :
			for transcript in dicoFeature[gene]["Transcript"] :
				#~ print str(gene) +"\t"+ str(transcript)
				for exon in dicoFeature[gene]["Transcript"][transcript]["Exon"]:
					chromosome = dicoFeature[gene]["Transcript"][transcript]["Exon"][exon]["Chromosome"]
					biotype = dicoFeature[gene]["Transcript"][transcript]["Exon"][exon]["Biotype"]
					start = dicoFeature[gene]["Transcript"][transcript]["Exon"][exon]["Start"]
					end = dicoFeature[gene]["Transcript"][transcript]["Exon"][exon]["End"]
					rank = dicoFeature[gene]["Transcript"][transcript]["Exon"][exon]["Rank"]
					strand = str(dicoFeature[gene]["Transcript"][transcript]["Exon"][exon]["Strand"])
					output.write(gene+"\t"+transcript+"\t"+chromosome+"\t"+biotype+"\t\t\t\t\t"+start+"\t"+end+"\t"+rank+"\t"+strand+"\n")
				if "5UTR" in dicoFeature[gene]["Transcript"][transcript] :
					chromosome = dicoFeature[gene]["Transcript"][transcript]["5UTR"]["Chromosome"]
					biotype = dicoFeature[gene]["Transcript"][transcript]["5UTR"]["Biotype"]
					start = dicoFeature[gene]["Transcript"][transcript]["5UTR"]["Start"]
					end = dicoFeature[gene]["Transcript"][transcript]["5UTR"]["End"]
					strand = str(dicoFeature[gene]["Transcript"][transcript]["5UTR"]["Strand"])
					output.write(gene+"\t"+transcript+"\t"+chromosome+"\t"+biotype+"\t"+start+"\t"+end+"\t\t\t\t\t\t"+strand+"\n")
				if "3UTR" in dicoFeature[gene]["Transcript"][transcript] :
					chromosome = dicoFeature[gene]["Transcript"][transcript]["3UTR"]["Chromosome"]
					biotype = dicoFeature[gene]["Transcript"][transcript]["3UTR"]["Biotype"]
					start = dicoFeature[gene]["Transcript"][transcript]["3UTR"]["Start"]
					end = dicoFeature[gene]["Transcript"][transcript]["3UTR"]["End"]
					strand = str(dicoFeature[gene]["Transcript"][transcript]["3UTR"]["Strand"])
					output.write(gene+"\t"+transcript+"\t"+chromosome+"\t"+biotype+"\t\t\t"+start+"\t"+end+"\t\t\t\t"+strand+"\n")
				
		output.close()
		
		#~ output = open("/home/anais/Documents/Data/"+sp+"/"+ini+"_gene_list.txt","w") # file opening for reading
		output = open("/home/anais/Documents/Data/Genomes/"+sp+"/"+ini+"_GeneID.txt","w") # file opening for reading
		for gene in dicoFeature :
			output.write(gene+"\n")
		output.close()
		print "Done"
	else :
		print "Echec for : "+sp

main()
