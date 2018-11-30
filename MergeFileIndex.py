#! usr/bin/python
import csv, math, numpy
import string
import subprocess
import Bio.Align.Applications
import sys
import Bio
import os, shutil
import pylab
import scipy.sparse as sparse
import re
import random
import matplotlib as mpl 
import matplotlib.pyplot as plt
import numpy as np
import argparse
import numpy as np 
import math 
from Bio import SeqIO
from scipy import stats
from pylab import *
import __future__
from math import sqrt
from numpy import mean, std

def AddG4(path, chromosome, specie, filename, suffixe, start_header):
	if (os.path.isfile(path+'/'+chromosome+'/'+specie+'_transcript_unspliced_'+chromosome+'_'+suffixe+'.txt') == True):
		input_file=path+'/'+chromosome+'/'+specie+'_transcript_unspliced_'+chromosome+'_'+suffixe+'.txt'
		#~ print(input_file)
		if (os.path.isfile(filename) == True): # if file out already exist
			output_file= open(filename,"a") ## file opening
		else: # if file out dont exist
			output_file= open(filename,"w") ## file opening
		
		inputfile= open(input_file,"r") ## file opening
		for line in inputfile:
			if (start_header not in line and start_header !=''): ## without header of files
				output_file.write(line)
			elif (start_header ==''):
				output_file.write(line)

		

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Count_G4')
	parser.add_argument ('-pi', '--pathInput', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/mouse')
	parser.add_argument ('-po', '--pathOutput', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/mouse')
	parser.add_argument ('-specie', '--specie', default = 'MM')
	return parser


def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	pathInput=arg.pathInput
	pathOutput=arg.pathOutput
	specie=arg.specie


	directory=pathOutput+'/All/'

	## Creation of directory for all chromosome of one specie
	#~ if os.path.isdir(directory):## if directory output exist 	
		#~ shutil.rmtree(directory) # delate if --> refresh
	#~ os.mkdir(directory) # create directory output

	filename=directory+specie+'_all_Index.txt'
	#lenghtFile=directory+specie+'_all_LenghtByGene.txt'
	#junctionFile=directory+specie+'_all_JunctionByGene.txt'
	#sectionFile=directory+specie+'_all_SectionByGene.txt'
	#transcriptFile=directory+specie+'_all_TranscriptPerG4.txt'
	#localisationFile=directory+specie+'_all_LocalisationPerTranscript.txt'
	Index=directory+specie+'_all_Index.txt'
	


	for path, dirs, files in os.walk(pathInput):
		for chromosome in dirs:
			if ('chr' in chromosome): # for each chromosome of this specie
				#if (os.path.isfile(path+'/'+chromosome+'/'+specie+'_'+chromosome+'_G4Uniqu.txt') == True):	
					#print path+'/'+chromosome+'/'+specie+'_'+chromosome+'_G4Uniqu.txt'
				#AddG4(pathInput, chromosome, specie, filename, 'G4Uniqu', 'InfoG4ByGene')
				#AddG4(pathInput, chromosome, specie, lenghtFile, 'LenghtByGene', 'gene_biotype')
				#AddG4(pathInput, chromosome, specie, junctionFile, 'JunctionByGene', 'gene_biotype')
				#AddG4(pathInput, chromosome, specie, sectionFile, 'SectionByGene', 'gene_biotype')
				#AddG4(pathInput, chromosome, specie, transcriptFile, 'TranscriptPerG4', 'Coordonnees')
				#AddG4(pathInput, chromosome, specie, localisationFile, 'LocalisationPerTranscript', 'InfoG4ByTranscript')
				#~ print(pathInput, chromosome, filename)
				AddG4(pathInput, chromosome, specie, filename, 'Index', '')

			
				 
main()

