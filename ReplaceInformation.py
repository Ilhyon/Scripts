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
import urllib 
# ~ import urllib2
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
			words=line.split(',')
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
	return exon_total

def SpecieSelction(sp):
	if sp == "APSH":
		return "Anaplasma_phagocytophilum_str_hz"
	elif sp == "AC":
		return "Anolis_carolinensis"
	elif sp == "AM":
		return "Apis_mellifera"
	elif sp == "AAV":
		return "Aquifex_aeolicus_vf5"
	elif sp == "AT":
		return "Arabidopsis_thaliana"
	elif sp == "AFD4":
		return "Archaeoglobus_fulgidus_dsm_4304"
	elif sp == "AN":
		return "Aspergillus_nidulans"
	elif sp == "BSSSS1":
		return "Bacillus_subtilis_subsp_subtilis_str_168"
	elif sp == "BBB":
		return "Borrelia_burgdorferi_b31"
	elif sp == "BAB1S99":
		return "Brucella_abortus_bv_1_str_9_941"
	elif sp == "CE":
		return "Caenorhabditis_elegans"
	elif sp == "CJSJN1A7":
		return "Campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819"
	elif sp == "CKCO":
		return "Candidatus_korarchaeum_cryptofilum_opf8"
	elif sp == "CSA":
		return "Cenarchaeum_symbiosum_a"
	elif sp == "CTDU3C":
		return "Chlamydia_trachomatis_d_uw_3_cx"
	elif sp == "CR":
		return "Chlamydomonas_reinhardtii"
	elif sp == "CAJ1F":
		return "Chloroflexus_aurantiacus_j_10_fl"
	elif sp == "CC":
		return "Chondrus_crispus"
	elif sp == "DR":
		return "Danio_rerio"
	elif sp == "DD":
		return "Dictyostelium_discoideum"
	elif sp == "DM":
		return "Drosophila_melanogaster"
	elif sp == "EH":
		return "Emiliania_huxleyi"
	elif sp == "EFV":
		return "Enterococcus_faecalis_v583"
	elif sp == "ECSK1SM":
		return "Escherichia_coli_str_k_12_substr_mg1655"
	elif sp == "FTSTSS":
		return "Francisella_tularensis_subsp_tularensis_schu_s4"
	elif sp == "GG":
		return "Gallus_gallus"
	elif sp == "GA":
		return "Gasterosteus_aculeatus"
	elif sp == "GSP":
		return "Geobacter_sulfurreducens_pca"
	elif sp == "HIRK":
		return "Haemophilus_influenzae_rd_kw20"
	elif sp == "HSR":
		return "Halobacterium_salinarum_r1"
	elif sp == "HS":
		return "Homo_sapiens"
	elif sp == "HBD5":
		return "Hyperthermus_butylicus_dsm_5456"
	elif sp == "LPSP":
		return "Legionella_pneumophila_str_paris"
	elif sp == "LM":
		return "Leishmania_major"
	elif sp == "MSA3":
		return "Methanobrevibacter_smithii_atcc_35061"
	elif sp == "MAC":
		return "Methanosarcina_acetivorans_c2a"
	elif sp == "MD":
		return "Monodelphis_domestica"
	elif sp == "MM":
		return "Mus_musculus"
	elif sp == "MTH":
		return "Mycobacterium_tuberculosis_h37rv"
	elif sp == "MPM":
		return "Mycoplasma_pneumoniae_m129"
	elif sp == "MXD1":
		return "Myxococcus_xanthus_dk_1622"
	elif sp == "NEKM":
		return "Nanoarchaeum_equitans_kin4_m"
	elif sp == "NMZ":
		return "Neisseria_meningitidis_z2491"
	elif sp == "NC":
		return "Neurospora_crassa"
	elif sp == "OA":
		return "Ornithorhynchus_anatinus"
	elif sp == "OS":
		return "Oryza_sativa"
	elif sp == "PT":
		return "Pan_troglodyte"
	elif sp == "PP":
		return "Physcomitrella_patens"
	elif sp == "PA":
		return "Pongo_abelii"
	elif sp == "PASI":
		return "Pyrobaculum_aerophilum_str_im2"
	elif sp == "PHO":
		return "Pyrococcus_horikoshii_ot3"
	elif sp == "SC":
		return "Saccharomyces_cerevisiae"
	elif sp == "SP":
		return "Schizosaccharomyces_pombe"
	elif sp == "SL":
		return "Solanum_lycopersicum"
	elif sp == "SASAN":
		return "Staphylococcus_aureus_subsp_aureus_n315"
	elif sp == "SPT":
		return "Streptococcus_pneumoniae_tigr4"
	elif sp == "SSP":
		return "Sulfolobus_solfataricus_p2"
	elif sp == "TAD1":
		return "Thermoplasma_acidophilum_dsm_1728"
	elif sp == "TTH":
		return "Thermus_thermophilus_hb8"
	elif sp == "VCOBETSN":
		return "Vibrio_cholerae_o1_biovar_el_tor_str_n16961"
	elif sp == "VV":
		return "Vitis_vinifera"
	elif sp == "WEODM":
		return "Wolbachia_endosymbiont_of_drosophila_melanogaster"
	elif sp == "YPBMS9":
		return 	"Yersinia_pestis_biovar_microtus_str_91001"

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
	nameSpecie = SpecieSelction(specie)
	ensembl = ("Homo_sapiens","Pan_troglodytes","Pongo_abelii","Mus_musculus","Monodelphis_domestica","Anolis_carolinensis","Ornithorhynchus_anatinus","Gallus_gallus","Danio_rerio","Gasterosteus_aculeatus")
	metazoa = ("Caenorhabditis_elegans","Drosophila_melanogaster","Apis_mellifera")
	plants = ("Oryza_sativa","Chlamydomonas_reinhardtii","Chondrus_crispus","Physcomitrella_patens","Solanum_lycopersicum","Vitis_vinifera","Arabidopsis_thaliana")
	fungi = ("Aspergillus_nidulans","Neurospora_crassa","Saccharomyces_cerevisiae","Schizosaccharomyces_pombe")
	protist = ("Dictyostelium_discoideum","Emiliania_huxleyi","Leishmania_major")
	bacteria = ("Mycoplasma_pneumoniae","Staphylococcus_aureus","Bacillus_subtilis","Enterococcus_faecalis","Streptococcus_pneumoniae","Chlamydia_trachomatis","Borrelia_burgdorferi","Mycobacterium_tuberculosis","Thermus_thermophilus","Geobacter_sulfurreducens","Campylobacter_jejuni","Wolbachia","Brucella_abortus","Escherichia_coli","Yersinia_pestis","Anaplasma_phagocytophilum","Aquifex_aeolicus","Chloroflexus_aurantiacus","Francisella_tularensis","Haemophilus_influenzae","Legionella_pneumophila","Myxococcus_xanthus","Neisseria_meningitidis","Vibrio_cholerae","Nanoarchaeum_equitans","Cenarchaeum_symbiosum","Pyrobaculum_aerophilum","Sulfolobus_solfataricus","Thermoplasma_acidophilum","Methanosarcina_acetivorans","Pyrococcus_horikoshii","Archaeoglobus_fulgidus","Candidatus_Korarchaeum","Halobacterium_salinarum","Hyperthermus_butylicus","Methanobrevibacter")
	archive_jul2016 = ("Gallus_gallus","Pan_trogolodytes")
	archive_feb2014 = ("Homo_sapiens","Mus_musculus")
	boo = False
	if specie in ensembl :
		url = 'http://useast.ensembl.org/'+nameSpecie+'/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;g='+geneID+';output=fasta;r='+chromosome+':'+str(start)+'-'+str(end)+';strand='+strand+';genomic=unmasked;_format=Text'
	elif specie in metazoa :
		url = 'https://metazoa.ensembl.org/'+nameSpecie+'/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;g='+geneID+';output=fasta;r='+chromosome+':'+str(start)+'-'+str(end)+';strand='+strand+';genomic=unmasked;_format=Text'
	elif specie in plants :
		url = 'https://plants.ensembl.org/'+nameSpecie+'/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;g='+geneID+';output=fasta;r='+chromosome+':'+str(start)+'-'+str(end)+';strand='+strand+';genomic=unmasked;_format=Text'
	elif specie in fungi :
		url = 'https://fungi.ensembl.org/'+nameSpecie+'/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;g='+geneID+';output=fasta;r='+chromosome+':'+str(start)+'-'+str(end)+';strand='+strand+';genomic=unmasked;_format=Text'
	elif specie in protist :
		url = 'https://protist.ensembl.org/'+nameSpecie+'/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;g='+geneID+';output=fasta;r='+chromosome+':'+str(start)+'-'+str(end)+';strand='+strand+';genomic=unmasked;_format=Text'
	elif specie in bacteria :
		url = 'https://bacteria.ensembl.org/'+nameSpecie+'/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;g='+geneID+';output=fasta;r='+chromosome+':'+str(start)+'-'+str(end)+';strand='+strand+';genomic=unmasked;_format=Text'
	else :
		print(specie)
	if specie in archive_jul2016 :
		url = 'http://jul2016.archive.ensembl.org/'+nameSpecie+'/Export/Output/Gene?db=core;flank3_display=0;flank5_display=0;g='+geneID+';output=fasta;r='+chromosome+':'+str(start)+'-'+str(end)+';strand='+strand+';genomic=unmasked;_format=Text'
		boo = True
	elif specie in archive_feb2014 :
		url = 'http://feb2014.archive.ensembl.org/'+nameSpecie+'/Export/Output/Gene?db=core;flank3_display=0;flank5_display=0;g='+geneID+';output=fasta;r='+chromosome+':'+str(start)+'-'+str(end)+';strand='+strand+';genomic=unmasked;_format=Text'
		boo = True
	if urllib2.urlopen(url) :
		response = urllib2.urlopen(url)
		fasta = response.read()
	else :
		print("Specie : "+specie+"\tgeneID : "+geneID)
		if boo :
			print(archive+"\n----------------")
	return fasta

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
	
def CreateSequence(directory, inputfilename, IntronPerGene, InfoPerGene, extension,specie):
	output1= open(directory+inputfilename.split(".")[0]+"_Sequence.txt","w") ## file opening
	for key, value in IntronPerGene.items():# for every intron by gene
		if (InfoPerGene.has_key(key) == True): 
			geneID=key
			chromosome=InfoPerGene.get(key).split("\t")[0]
			biotype=InfoPerGene.get(key).split("\t")[1]
			strand=InfoPerGene.get(key).split("\t")[2]
			for value_intron in value:
				start_intron=value_intron[0]
				end_intron=value_intron[1]
				if (strand == str(1)): ## if gene positif
					start_sup=int(start_intron)-1-extension
					end_sup=int(start_intron)-1
					start_inf=int(end_intron)+1
					end_inf=int(end_intron)+1+extension
					fasta_amont=Fasta(geneID,chromosome,start_sup,end_sup,strand,specie)
					sequence_amont=Sequence(fasta_amont, strand)
					#if (sequence_amont):
						#print' yes'
					#else:
						#print ' no'
					fasta_aval=Fasta(geneID,chromosome,start_inf,end_inf,strand,specie)
					sequence_aval=Sequence(fasta_aval, strand)
			
				else :
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
				#print sequence
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
			print("error exon "+geneID+" "+transcriptID)
		if (IntronPerTranscript.has_key(transcriptID) == True): 
				intron=IntronPerTranscript.get(transcriptID)
				for couple in intron:
					intronList=intronList+'-'.join(couple)+";"
		else:
			print("error intron "+geneID+" "+transcriptID)
		
		output5.write(transcriptID+"|"+geneID+"|"+chromosome+"|"+strand+"|"+biotype+"|"+exonList[:-1]+"|"+intronList[:-1]+"|"+start5+"|"+end5+"|"+start3+"|"+end3+"\n")
	
		
def build_arg_parser():
	GITDIR=os.getcwd()+'/../'
	parser = argparse.ArgumentParser(description = 'ReplaceInformationBiomart')
	parser.add_argument ('-p', '--path', default = GITDIR+'Data/Genomes/')
	parser.add_argument ('-sp', '--specie', default = 'HS')
	parser.add_argument ('-ext', '--extension', default = 100)
	return parser

def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	directory=arg.path
	sp=arg.specie
	extension=arg.extension
	
	words = sp.split("_")
	letters = [word[0] for word in words]
	specie = "".join(letters)
	specie = ini.upper()
	
	inputfilename = specie+'_transcript_unspliced.txt'
	
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
	
###################################################################################################################################### CREATE FILES
	##Create file Sequence
	CreateSequence(directory, inputfilename,IntronPerGene, InfoPerGene, extension,specie)
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

	CreateIndex(directory, inputfilename,InformationPerGeneAndTranscript,ExonPerTranscript, IntronPerTranscript)

	
				
		
			
	
main()
