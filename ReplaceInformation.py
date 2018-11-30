#! usr/bin/python

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
import numpy as np
import argparse
from Bio import SeqIO
from scipy import stats
import copy
import urllib 
import urllib2
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
		if (re.search('ENS', line)): 
			words=line.split('\t')
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
	exonTotalOrdered=[]
	while i < len(rank):
			exon=[]
			if (strand == str(1)): ## if gene positif
				exon.append(str(start_exon[i]))
				exon.append(str(end_exon[i]))
				i=i+1
			else:  ## if gene negatif
				exon.append(str(end_exon[i]))
				exon.append(str(start_exon[i]))
				i=i+1
			exon_total.append(exon) ## create array of array with start and end of each exon 
	for i in range(len(rank)):
		exonTotalOrdered.append(0)
	for i in range(len(rank)):
		numberExon=rank[i]
		exonTotalOrdered[int(numberExon)-1]=exon_total[i]
	return exonTotalOrdered

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

def CreateEndIntron(geneID,exon_total, rank, strand):
	end_intron=[]
	
	for i in exon_total:
		
		if ( int(i[0]) != int(exon_total[-len(rank)][0]) ): ## EXCEPT start firt exon
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

def Fasta(geneID,chromosome,start,end, strand,specie):
	nameSpecie=''
	if (specie == 'HS'): # human
		nameSpecie='Homo_sapiens'
	elif (specie == 'MM'):# mouse
		nameSpecie='Mus_musculus'
	url = 'http://useast.ensembl.org/'+nameSpecie+'/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;g='+geneID+';output=fasta;r='+chromosome+':'+str(start)+'-'+str(end)+';strand='+strand+';genomic=unmasked;_format=Text'
	response = urllib2.urlopen(url)
	fasta = response.read()
	return fasta

def Sequence(fasta, strand):
	sequence=""
	for i in fasta.split("\r\n"):
		if not (re.search('^>', i) or re.search('\n', i)):
			sequence=sequence+i
	return sequence
	
def CreateSequence(directory, inputfilename, IntronPerGene, InfoPerGene, extension,specie):
	output1= open(directory+inputfilename.split(".")[0]+"_Sequence.txt","w") ## file opening
	for key, value in IntronPerGene.items():# for every intron by gene
		if (InfoPerGene.has_key(key) == True): 
			geneID=key
			chromosome=InfoPerGene.get(key).split("\t")[0]
			biotype=InfoPerGene.get(key).split("\t")[1]
			strand=InfoPerGene.get(key).split("\t")[2]
			for value_intron in value: # for each intron of this gene
				start_intron=value_intron[0]
				end_intron=value_intron[1]
				if (strand == str(1)): ## if gene positif
					start_sup=int(start_intron)-1-extension
					end_sup=int(start_intron)-1
					start_inf=int(end_intron)+1
					end_inf=int(end_intron)+1+extension
					fasta_amont=Fasta(geneID,chromosome,start_sup,end_sup,strand,specie)
					sequence_amont=Sequence(fasta_amont, strand)
					fasta_aval=Fasta(geneID,chromosome,start_inf,end_inf,strand,specie)
					sequence_aval=Sequence(fasta_aval, strand)
			
				else : ## if gene negatif
					start_sup=int(start_intron)+1+extension
					end_sup=int(start_intron)+1
					start_inf=int(end_intron)-1
					end_inf=int(end_intron)-1-extension
					fasta_amont=Fasta(geneID,chromosome,end_sup,start_sup,strand,specie)
					#print fasta_amont
					sequence_amont=Sequence(fasta_amont, strand)
					fasta_aval=Fasta(geneID,chromosome,end_inf,start_inf,strand,specie)
					sequence_aval=Sequence(fasta_aval, strand)
				sequence=sequence_amont+sequence_aval
				output1.write(">"+geneID+"|"+"-".join(value_intron)+"\n"+str(sequence)+"\n")

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
			print "error exon", geneID, transcriptID
		if (IntronPerTranscript.has_key(transcriptID) == True): 
				intron=IntronPerTranscript.get(transcriptID)
				for couple in intron:
					intronList=intronList+'-'.join(couple)+";"
		else:
			print "error intron", geneID, transcriptID
		
		output5.write(transcriptID+"|"+geneID+"|"+chromosome+"|"+strand+"|"+biotype+"|"+exonList[:-1]+"|"+intronList[:-1]+"|"+start5+"|"+end5+"|"+start3+"|"+end3+"\n")
	
		
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'ReplaceInformationBiomart')
	parser.add_argument ('-p', '--path', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/mouse')
	parser.add_argument ('-chr', '--chromosome', default = 'X')
	parser.add_argument ('-specie', '--specie', default = 'MM')
	parser.add_argument ('-ext', '--extension', default = 100)
	return parser



def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path
	chromosome=arg.chromosome
	specie=arg.specie
	extension=arg.extension
	
	directory=path+'/chr'+chromosome+'/'
	inputfilename = specie+'_transcript_unspliced_chr'+chromosome+'.txt'
	
	
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

		if (re.compile('[1-9, X, Y]' ).search(chromosome) and not re.compile('[A-W, Z]' ).search(chromosome)):
				InfoPerGene[geneID] = chromosome+'\t'+biotype+'\t'+strand ### ajout de 5' et 3' ?
											### if biotype gene = coding mais absence de 5' et 3' --> biotype transcript = non coding
		exon_total=ExonTotal(rank, strand, start_exon, end_exon)
		start_intron=CreateStartIntron(exon_total, rank, strand)
		end_intron=CreateEndIntron(geneID,exon_total, rank, strand)
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
	
	
	
###################################################################################################################################### CREATE FILES
	##Create file Sequence
	
	#CreateSequence(directory, inputfilename,IntronPerGene, InfoPerGene, extension,specie)
	

	CreateIndex(directory, inputfilename,InformationPerGeneAndTranscript,ExonPerTranscript, IntronPerTranscript)

	
	
				
		
			
	
main()
