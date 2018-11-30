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

def AddInfo(listFiles, pathInput, specie, filename, suffixe, start_header):
	for fileInDirectory in listFiles:
		if (suffixe in fileInDirectory and 'Unique' not in fileInDirectory and '.fas' not in fileInDirectory):
			#~ print fileInDirectory
			if (os.path.isfile(filename) == True): # if file out already exist
				output_file= open(filename,"a") ## file opening
			else: # if file out dont exist
				output_file= open(filename,"w") ## file opening
			inputfile= open(pathInput+'/'+fileInDirectory,"r") ## file opening
			for line in inputfile:
				if (start_header not in line and start_header !=''): ## without header of files
					output_file.write(line)
				elif (start_header ==''):
					output_file.write(line)

		

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Count_G4')
	parser.add_argument ('-pi', '--pathInput', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/mouse/G4Screener')
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
	directoryTranscriptType = pathInput+'/../transcriptType'

	## Creation of directory for all chromosome of one specie
	if os.path.isdir(directory):## if directory output exist 	
		shutil.rmtree(directory) # delate if --> refresh
	os.mkdir(directory) # create directory output

	TranscriptPerG4=directory+specie+'_all_TranscriptPerG4.txt'
	G4InGenome=directory+specie+'_all_G4InGenome.txt'
	G4InTranscript=directory+specie+'_all_G4InTranscript.txt'
	TranscriptType=directory+specie+'_all_TranscriptType.txt'
	
	listFiles = os.listdir(pathInput)
	listFilesT = os.listdir(directoryTranscriptType)
	AddInfo(listFiles, pathInput, specie, TranscriptPerG4, 'TranscriptPerG4', 'Transcript(s)')
	AddInfo(listFiles, pathInput, specie, G4InGenome, 'G4InGenome', 'InfoG4ByGene')
	AddInfo(listFiles, pathInput, specie, G4InTranscript, 'G4InTranscript', 'InfoG4ByTranscript')
	AddInfo(listFilesT, directoryTranscriptType, specie, TranscriptType, 'transcriptType', '')

				 
main()

