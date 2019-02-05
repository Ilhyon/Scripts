import csv, math, numpy
import string
import subprocess
import Bio.Align.Applications
import sys
import Bio
import os
import scipy.sparse as sparse
import re
import random
# ~ import matplotlib.pyplot as plt
import numpy as np
import argparse
from Bio import SeqIO
from scipy import stats
import copy
from Bio import SeqIO

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import gzip

from cogent.db.ensembl import HostAccount
from cogent.db.ensembl import Species
#from pyGeno.Genome import *

def OrderInformationBiomart(directory,inputfilename): ### create one line per transcript, start and end no take information of strand here 
	old_transcript=""
	InformationPerGeneAndTranscript={}
	inputfile= open(directory+inputfilename,"r") 
	for line in inputfile:
		if not (re.search('Gene stable ID', line)): 
			words=line.split("\t")
			# ~ print words
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
	return exon_total

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
			
def CreateStartIntron(exon_total, rank, strand):
	start_intron=[]
	for i in exon_total:
		if ( int(i[1]) != int(exon_total[-1][1]) ): ## EXCEPT end last exon
			if (strand == str(1)):
				debut_intron=int(i[1])+1 ## end exon +1 == > debut intron 
			else:
				debut_intron=int(i[1])-1 ## end exon +1 == > debut intron 
			start_intron.append(debut_intron)
	return start_intron ## create array  with start of each intron 

def CreateEndIntron(exon_total, rank, strand):
	end_intron=[]
	for i in exon_total:
		if ( int(i[0]) != int(exon_total[-len(rank)][0]) ): ## EXCEPT start first exon
			if (strand == str(1)):
				fin_intron=int(i[0])-1 ## end exon +1 == > debut intron 
			else:
				fin_intron=int(i[0])+1
			end_intron.append(fin_intron)
	return end_intron ## create array  with end of each intron 

def IntronTotal(start_intron, end_intron, strand):
	i=0
	intron_total=[]
	while i < len(start_intron): ## for number of intron
			intron=[]
			intron.append(str(start_intron[i]))
			intron.append(str(end_intron[i]))
			intron_total.append(intron)
			i=i+1
	return intron_total ## create array of array with start and end of each intron 

def AddIntron(intron_total, Intron):
	for i in intron_total:
		if i not in Intron:
			Intron.append(i) 
	return Intron

def AddTranscriptPerIntron(intron_total, Dico, transcriptID):
	for i in intron_total:
		if (Dico.has_key(i[0]+'-'+i[1]) == False): ## if this intron is not in dico 
			Dico[i[0]+'-'+i[1]] = transcriptID
		else:
			if bool(Dico): ## dico already with info
				last_transcript= Dico.get(i[0]+'-'+i[1])
				Dico[i[0]+'-'+i[1]] = last_transcript +"-"+transcriptID
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
	sequence = dicoChromosome[chromosome][(int(start)-1):(int(end)-1+1)]
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

def Sequence(fasta, strand):
	sequence=""
	for i in fasta.split("\r\n"):
		if not (re.search('^>', i) or re.search('\n', i)):
			#if (strand == str(1)):
			sequence=sequence+i
			#else:
				#dna = Seq(i, generic_dna)
				#sequence_reverse= dna.reverse_complement()
				#sequence=sequence_reverse+sequence
				#sequence=i+sequence
				
	return sequence
	
def CreateSequence(directory, inputfilename, IntronPerGene, InfoPerGene, extension,specie,nameSpecie):
	output1= open(directory+inputfilename.split(".")[0]+"_Sequence.txt","w") ## file opening
	for key, value in IntronPerGene.items():# for every intron by gene
		if key in InfoPerGene :
			geneID=key
			chromosome=InfoPerGene.get(key).split("\t")[0]
			biotype=InfoPerGene.get(key).split("\t")[1]
			strand=InfoPerGene.get(key).split("\t")[2]
			DicoSequence = ImportFasta(nameSpecie)
			for value_intron in value:
				start_intron=int(value_intron[0])
				end_intron=int(value_intron[1])
				if (strand == str(1)): ## if gene positif
					start_sup=start_intron-1-extension
					end_sup=start_intron-1
					start_inf=end_intron+1
					end_inf=end_intron+1+extension
					fasta_amont=Fasta(chromosome,start_sup,end_sup,strand,nameSpecie,DicoSequence)
					sequence_amont=Sequence(fasta_amont, strand)
					fasta_aval=Fasta(chromosome,start_inf,end_inf,strand,nameSpecie,DicoSequence)
					sequence_aval=Sequence(fasta_aval, strand)
			
				else :
					start_sup=start_intron+1+extension
					end_sup=start_intron+1
					start_inf=end_intron-1
					end_inf=end_intron-1-extension
					fasta_amont=Fasta(chromosome,end_sup,start_sup,strand,nameSpecie,DicoSequence)
					#print fasta_amont
					sequence_amont=Sequence(fasta_amont, strand)
					fasta_aval=Fasta(chromosome,end_inf,start_inf,strand,nameSpecie,DicoSequence)
					sequence_aval=Sequence(fasta_aval, strand)
				sequence=sequence_amont+sequence_aval
				#print sequence
				if geneID == "Cj1325":
					print value_intron
					print fasta_amont
				output1.write(">"+geneID+"|"+"-".join(value_intron)+"\n"+str(sequence)+"\n")
	output1.close()
	print "\t\t Junction done."

def CreateIndex(directory, inputfilename,InformationPerGeneAndTranscript,ExonPerTranscript, IntronPerTranscript):
	output5= open(directory+inputfilename.split(".")[0]+"_Index.txt","w")
	for key, values in InformationPerGeneAndTranscript.items():
		exonList=""
		intronList=""
		geneID=key.split("-")[0]
		transcriptID=key.split("-")[1]
		chromosome=values.split("|")[0]
		biotype=values.split("|")[1]
		start5=values.split("|")[2]
		end5=values.split("|")[3]
		start3=values.split("|")[4]
		end3=values.split("|")[5]
		start_exon=values.split("|")[6].split(";")
		strand=values.split("|")[9]
		if (ExonPerTranscript.has_key(transcriptID) == True): 
				exon=ExonPerTranscript.get(transcriptID)
				for couple in exon:
					exonList=exonList+'-'.join(couple)+";"
		else:
			print("error exon "+geneID+" "+transcriptID)
		if (IntronPerTranscript.has_key(transcriptID) == True): 
				intron=IntronPerTranscript.get(transcriptID)
				for couple in intron:
					intronList=intronList+'-'.join(couple)+";"
		else:
			print("error intron "+geneID+" "+transcriptID)
		output5.write(transcriptID+"|"+geneID+"|"+chromosome+"|"+strand+"|"+biotype+"|"+exonList[:-1]+"|"+intronList[:-1]+"|"+start5+"|"+end5+"|"+start3+"|"+end3+"\n")
	output5.close()
	print "\t\t Index done."
		
def build_arg_parser():
	GITDIR=os.getcwd()+'/../'
	parser = argparse.ArgumentParser(description = 'ReplaceInformationBiomart')
	parser.add_argument ('-p', '--path', default = "/home/anais/Documents/Data/Genomes/")
	parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
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
	
	
	for key, values in InformationPerGeneAndTranscript.items():
		geneID=key.split("-")[0]
		transcriptID=key.split("-")[1]
		chromosome=values.split("|")[0]
		biotype=values.split("|")[1]
		start5=values.split("|")[2]
		end5=values.split("|")[3]
		start3=values.split("|")[4]
		end3=values.split("|")[5]
		start_exon=values.split("|")[6].split(";")
		end_exon=values.split("|")[7].split(";")
		rank=values.split("|")[8].split(";")
		strand=values.split("|")[9]

		# ~ if (re.compile('[1-9, X, Y]' ).search(chromosome) and not re.compile('[A-W, Z]' ).search(chromosome)):
		InfoPerGene[geneID] = chromosome+'\t'+biotype+'\t'+strand ### ajout de 5' et 3' ?
											### if biotype gene = coding mais absence de 5' et 3' --> biotype transcript = non coding
		exon_total=ExonTotal(rank, strand, start_exon, end_exon)
		exon_total = sortExon(exon_total)
		start_intron=CreateStartIntron(exon_total, rank, strand)
		end_intron=CreateEndIntron(exon_total, rank, strand)
		intron_total = IntronTotal(start_intron, end_intron, strand)
		ExonPerTranscript[transcriptID] = exon_total
		IntronPerTranscript[transcriptID] = intron_total
		if (old_gene != geneID): ## if new gene	
			Intron= []  ## create new Intron
			Intron=AddIntron(intron_total, Intron)	## copy information transcript #1
			old_gene=geneID
		else: ## if other transcript
			Intron=AddIntron(intron_total, Intron)	 ## copy information other transcripts
			old_gene=geneID
		
		IntronPerGene[geneID] = Intron ## add information in dico (modification of value if key already exist)
		TranscriptPerIntron=AddTranscriptPerIntron(intron_total, TranscriptPerIntron, transcriptID)
		# ~ if transcriptID == "CAL35439":
			# ~ print IntronPerGene
			
###################################################################################################################################### CREATE FILES
	##Create file Sequence
	CreateSequence(directory,inputfilename,IntronPerGene,InfoPerGene, extension,specie,nameSpecie)
	'''
	##Create file IntronPerGene		
	output2= open(directory+inputfilename.split(".")[0]+"_IntronPerGene.txt","w") 	
	for key, value in IntronPerGene.items():
		string=""
		for i in value:
			string=string+"-".join(i)+"|"
		output2.write(key+"\t"+string+"\n")

	##Create file TranscriptPerIntron
	output3= open(directory+inputfilename.split(".")[0]+"_TranscriptPerIntron.txt","w") 	
	for key, value in TranscriptPerIntron.items():
		output3.write(key+"\t"+value+"\n")

	##Create file InfoPerGene
	output4= open(directory+inputfilename.split(".")[0]+"_InfoPerGene.txt","w") 	
	for key, value in InfoPerGene.items():
		output4.write(key+"\t"+value+"\n")
	'''

	CreateIndex(directory,inputfilename,InformationPerGeneAndTranscript,ExonPerTranscript,IntronPerTranscript)

	
				
		
			
	
main()
