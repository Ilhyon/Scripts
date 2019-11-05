#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import sys
import Bio
import os
import re
from pprint import pprint
import argparse

def OrderInformationBiomart(directory,inputfilename): ### create one line per transcript, start and end no take information of strand here
	"""
		Return a dictionary of strings : {geneID-trID : chr|biotype|start5|end5|start3|end3|allStartExon|allEndExon|allRankExon|strand}
	"""
	old_transcript=""
	InformationPerGeneAndTranscript={}
	inputfile= open(directory+inputfilename,"r")
	for line in inputfile:
		if not (re.search('Gene stable ID', line)):
			words=line.split("\t")
			geneID=words[0].rstrip() # gene Id
			transcriptID=words[1].rstrip() # transcropt Id
			chromosome=words[2].rstrip()
			biotype=words[3].rstrip() ## biotype of the GENE
			start_5= words[4].rstrip()
			end_5= words[5].rstrip()
			start_3= words[6].rstrip()
			end_3= words[7].rstrip()
			start_exon= words[8].rstrip() # start for ONE exon
			end_exon= words[9].rstrip()
			rank= words[10].rstrip()
			strand= words[11].rstrip()
			# ~ if strand == "-1" and int(start_exon) > int(end_exon):
				# ~ start_exon= words[9].rstrip() # start for ONE exon
				# ~ end_exon= words[8].rstrip()
			if (old_transcript != transcriptID ): ## if new transcript
				##maj value
				var_rank= ""
				var_startExon= ""
				var_endExon=""
				start5=""
				end5=""
				start3=""
				end3=""
				# assignation
				var_rank=rank
				var_startExon=start_exon
				var_endExon=end_exon
				start5=start_5
				end5=end_5
				start3=start_3
				end3=end_3
				old_transcript = transcriptID
			else: ## add information for transcript with several exons
				var_rank=var_rank+";"+rank
				var_startExon=var_startExon+";"+start_exon
				var_endExon=var_endExon+";"+end_exon
				start3=start_3
				end3=end_3
			InformationPerGeneAndTranscript[geneID+"-"+transcriptID]=chromosome+"|"+biotype+"|"+start5+"|"+end5+"|"+start3+"|"+end3+"|"+var_startExon +"|"+var_endExon+"|"+var_rank+"|"+strand ## la
	return InformationPerGeneAndTranscript

def ExonTotal(rank, strand, start_exon, end_exon):
	"""
		rank, start_exon and end_exon -> list that contain all features for one transcript
		Return a list of exon list : [[Start exon 1 - end exon],...]
	"""
	i=0
	exon_total=[]
	while i < len(rank):
			exon=[]
			if (strand == str(1)): ## forward strand
				exon.append(str(start_exon[i]))
				exon.append(str(end_exon[i]))
				i=i+1
			else:  ## reverse strand
				exon.append(str(end_exon[i]))
				exon.append(str(start_exon[i]))
				i=i+1
			exon_total.append(exon) ## create array of array with start and end of each exon

	return sortExon(exon_total)

def sortExon(exon_total):
	exonSorted = []
	for exon in exon_total:
		if not exonSorted :
			exonSorted = [exon]
		elif exon[0] < exonSorted[0][0]: # the new exon start is the first exon
			exonSorted.insert(0, exon) # then we change the first exon
		elif exon[0] > exonSorted[-1][0]: # the new exon start is the last exon
			exonSorted.insert(-1, exon)
		else:
			for i in exon_total[1::-2]:
				if exon[0] < i[0]:
					exonSorted.insert(exon_total.index(i),exon)
				else:
					exonSorted.insert(exon_total.index(i)+1,exon)
	return exonSorted

def CreateStartIntron(exon_total, strand):
	"""
		Return a list with all the starts of intron for one transcript
	"""
	start_intron=[]
	lastEndExon = int(exon_total[-1][1])
	for coupleSE in exon_total:
		endExon = int(coupleSE[1])
		if endExon != lastEndExon:
			if strand == "1":
				debut_intron = endExon + 1 ## end exon +1 == > debut intron
			else:
				debut_intron= endExon - 1 ## end exon -1 == > debut intron
			start_intron.append(debut_intron)
	return start_intron

def CreateEndIntron(exon_total, strand):
	"""
		Return a list with all the ends of intron for one transcript
	"""
	end_intron=[]
	firstStartExon = int(exon_total[0][0])
	for coupleSE in exon_total:
		startExon = int(coupleSE[0])
		if startExon != firstStartExon : ## EXCEPT start first exon
			if strand == "1":
				fin_intron = startExon - 1 ## end exon +1 == > debut intron
			else:
				fin_intron = startExon + 1
			end_intron.append(fin_intron)
	return end_intron ## create array  with end of each intron

def IntronTotal(exons_total, strand):
	"""
		With the list of exon for one trancript, it generate a list of
		all introns : [[Start intron 1 - end intron],...].
	"""
	i=0
	introns_total=[]
	starts_intron = CreateStartIntron(exons_total, strand)
	ends_intron = CreateEndIntron(exons_total, strand)
	while i < len(starts_intron): ## for number of intron
			intron=[]
			intron.append(starts_intron[i])
			intron.append(ends_intron[i])
			introns_total.append(intron)
			i=i+1
	return introns_total ## create array of array with start and end of each intron

def AddIntron(introns_total, Intron):
	for coupleCoord in introns_total:
		if coupleCoord not in Intron:
			Intron.append(coupleCoord)
	return Intron

def AddTranscriptPerIntron(introns_total, Dico, transcriptID):
	"""
		Return a dictionary containing for each intron the transcripts
		where they are
		{StartIntron-EndIntron : trId1-trId2...}
	"""
	for intron in introns_total:
		intronCoord = str(intron[0])+'-'+str(intron[1])
		if intronCoord not in Dico: # this intron is not in dico
			Dico[intronCoord] = transcriptID
		else : # the intron is in the dico
			if Dico: ## dico already with info
				last_transcript = Dico[intronCoord]
				Dico[intronCoord] = last_transcript +"-"+transcriptID
	# ~ pprint(Dico)
	return Dico

def ImportFasta(nameSpecie):
	directory = "/home/anais/Documents/Data/Genomes/"+nameSpecie+"/Fasta/"
	dicoChromosome = {}
	listFile = os.listdir(directory)
	for filename in listFile :
		with open(directory+filename) as f: # file opening
			content = f.read()
			l = content.split('\n')
			if l[0].startswith('>'):
				header = l[0]
				chromosome = header.split(' ')[0][1:]
				sequence = "".join(l[1:])
				dicoChromosome.update({chromosome : sequence})
	return dicoChromosome

def Fasta(chromosome,start,end,strand,nameSpecie,dico):
	if chromosome in dico :
		sequence = dico[chromosome][(int(start)-1):(int(end)-1+1)]
	else :
		sequence = ""
		# ~ print nameSpecie+" have a probleme with the chromosome : "+chromosome
	if strand == "-1" :
		reverse = ""
		for n in sequence:
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
		sequence = reverse[::-1]
	return sequence

def CreateSequence(directory, inputfilename, IntronPerGene, InfoPerGene, extension,specie,nameSpecie):
	output= open(directory+inputfilename.split(".")[0]+"_Sequence.txt","w") ## file opening
	DicoSequence = ImportFasta(nameSpecie)
	print "\t\t Succeed to import the chromosomes sequences."
	for geneID in IntronPerGene: # for every gene
		chromosome = InfoPerGene[geneID].split("\t")[0]
		biotype = InfoPerGene[geneID].split("\t")[1]
		strand = InfoPerGene[geneID].split("\t")[2]
		for value_intron in IntronPerGene[geneID]: # browse all intron of the gene
			start_intron = value_intron[0]
			end_intron = value_intron[1]
			if strand == "1" : ## if gene positif
				start_sup = start_intron - 1 - extension
				end_sup = start_intron -1
				start_inf = end_intron +1
				end_inf = end_intron + 1 + extension
				sequence_amont = Fasta(chromosome,start_sup,end_sup,strand,nameSpecie,DicoSequence)
				sequence_aval = Fasta(chromosome,start_inf,end_inf,strand,nameSpecie,DicoSequence)

			else :
				start_sup = start_intron + 1 + extension
				end_sup = start_intron +1
				start_inf = end_intron -1
				end_inf = end_intron - 1 - extension
				sequence_amont = Fasta(chromosome,end_sup,start_sup,strand,nameSpecie,DicoSequence)
				sequence_aval = Fasta(chromosome,end_inf,start_inf,strand,nameSpecie,DicoSequence)
			sequence=sequence_amont+sequence_aval
			output.write(">"+geneID+"|"+"-".join(map(str,value_intron))+"\n"+str(sequence)+"\n")
	output.close()
	print "\t\t Junction done."

def build_arg_parser():
	GITDIR=os.getcwd()+'/../'
	parser = argparse.ArgumentParser(description = 'ReplaceInformationBiomart')
	parser.add_argument ('-p', '--path', default = "/home/anais/Documents/Data/Genomes/")
	parser.add_argument ('-sp', '--specie', default = 'pyrococcus_horikoshii_ot3')
	parser.add_argument ('-ext', '--extension', default = 100)
	return parser

def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	directory=arg.path
	nameSpecie=arg.specie
	extension=arg.extension

	words = nameSpecie.split("_")
	letters = [word[0] for word in words]
	specie = "".join(letters)
	specie = specie.upper()

	print "Parser started for : "+nameSpecie+"\t"+specie
	inputfilename = nameSpecie+"/"+specie+'_transcript_unspliced.txt'

	old_gene=""

	Intron=[]
	IntronPerGene = {}
	ExonPerTranscript = {}
	IntronPerTranscript = {}
	TranscriptPerIntron = {}
	InfoPerGene={}

	InformationPerGeneAndTranscript=OrderInformationBiomart(directory, inputfilename)

	output= open(directory+inputfilename.split(".")[0]+"_Index.txt","w")
	for idFeature in InformationPerGeneAndTranscript :
		geneID=idFeature.split("-")[0]
		transcriptID=idFeature.split("-")[1]
		info = InformationPerGeneAndTranscript[idFeature]
		chromosome=info.split("|")[0]
		if nameSpecie == "gasterosteus_aculeatus" and chromosome != "MT":
			# their chromosome name are like "groupI" so I only retrieve
			# the number of the chromosome
			chromosome = chromosome.split("p")[1]
		biotype=info.split("|")[1]
		start5=info.split("|")[2]
		end5=info.split("|")[3]
		start3=info.split("|")[4]
		end3=info.split("|")[5]
		starts_exon=info.split("|")[6].split(";")
		ends_exon=info.split("|")[7].split(";")
		ranks=info.split("|")[8].split(";")
		strand=info.split("|")[9]
		if len(chromosome) <= 3 : # it avoid to have contig with name like JH375512.1
			InfoPerGene[geneID] = chromosome+'\t'+biotype+'\t'+strand ### ajout de 5' et 3' ?
				### if biotype gene = coding mais absence de 5' et 3' --> biotype transcript = non coding
			exons_total = ExonTotal(ranks, strand, starts_exon, ends_exon)
			introns_total = IntronTotal(exons_total, strand)
			ExonPerTranscript[transcriptID] = exons_total
			IntronPerTranscript[transcriptID] = introns_total
			if (old_gene != geneID): # if new gene
				Intron = introns_total	## create new Intron and copy information transcript #1
				old_gene=geneID
			else: ## if other transcript
				Intron=AddIntron(introns_total, Intron)	 ## copy information other transcripts
				old_gene=geneID

			IntronPerGene[geneID] = Intron ## add information in dico (modification of value if key already exist)
			TranscriptPerIntron=AddTranscriptPerIntron(introns_total, TranscriptPerIntron, transcriptID)

			exon=ExonPerTranscript[transcriptID]
			exonList=""
			intronList=""
			for couple in exon:
				exonList=exonList+'-'.join(couple)+";"
			intron=IntronPerTranscript.get(transcriptID)
			for couple in intron:
				intronList=intronList+'-'.join(map(str,couple))+";"

			output.write(transcriptID+"|"+geneID+"|"+chromosome+"|"+strand+"|"+biotype+"|"+exonList[:-1]+"|"+intronList[:-1]+"|"+start5+"|"+end5+"|"+start3+"|"+end3+"\n")
	output.close()
	print "\t\t Index done."
	CreateSequence(directory,inputfilename,IntronPerGene,InfoPerGene, extension,specie,nameSpecie)

main()
