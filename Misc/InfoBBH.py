#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse

#----------------------------------------------------------------------#
def ImportG4sSpecie(filename):
	""" Create a list containing all the gene with at least one rG4s
	Parameters
	    ----------
	     filename : string, name of the file containing all the gene
	     with a rG4 in it
	    Returns
	    -------
	     listrG4s : list with the gene
	"""
	inputFile = open(filename,"r")
	listrG4s = []
	for line in inputFile:
		listrG4s.append(line)
	inputFile.close()
	return(listrG4s)
#----------------------------------------------------------------------#
def ImportBBH(filename):
	""" Create two list containing all the genes that are BBH between our
	two specie, so there is one list per specie
	Parameters
	    ----------
	     filename : string, name of the file containing all the gene
	     of both specie as BBH
	    Returns
	    -------
	     listBBHSpecie1 : list of string, list with the gene of the first specie
	     listBBHSpecie2 : list of string, list with the gene of the second specie
	"""
	inputFile = open(filename,"r")
	listBBHSpecie1 = []
	listBBHSpecie2 = []
	for line in inputFile:
		words=line.split('\t')
		listBBHSpecie1.append(words[0])
		listBBHSpecie2.append(words[1])
	inputFile.close()
	return(listBBHSpecie1,listBBHSpecie2)
#----------------------------------------------------------------------#
def GetBBHandrG4s(listBBH, listrG4s):
	""" Create a list of the genes that contains a rG4s and that is also
	a BBH with the other specie
	Parameters
	    ----------
	     listBBH : see ImportBBH
	     listrG4s : see ImportG4sSpecie
	    Returns
	    -------
	     listrG4sandBBH : list of string containing only the genes that have
	     a rG4s and that are a BBH
	"""
	listrG4sandBBH = []
	for gene in listBBH: # browse of the list of BBH and not rG4s
	#because it is the shortest and the result will be the same
		if gene in listrG4s:
			listrG4sandBBH.append(gene)
	return(listrG4sandBBH)
#----------------------------------------------------------------------#
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'InfoBBH')
	parser.add_argument ('-p', '--path', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/Homologie/BBH/')
	parser.add_argument ('-sp1', '--specie1', default = 'MM')
	parser.add_argument ('-sp2', '--specie2', default = 'HS')
	return parser
#----------------------------------------------------------------------#
def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path	# directory which contain all the file with BBH and rG4s
	specie1=arg.specie1	# first specie to analyse
	specie2=arg.specie2	# first specie to analyse
	
	# specie one
	filename = path+specie1+'_gene_rG4s.txt'
	specie1rG4s = ImportG4sSpecie(filename)
	# specie two
	filename = path+specie2+'_gene_rG4s.txt'
	specie2rG4s = ImportG4sSpecie(filename)
	
	# BBH
	filename = path+'BBH.txt'
	BBHspecie1, BBHspecie2 = ImportBBH(filename)
	
	#gene common
	BBHG4sSpecie1 = GetBBHandrG4s(BBHspecie1, specie1rG4s)
	BBHG4sSpecie2 = GetBBHandrG4s(BBHspecie2, specie2rG4s)
	
	# write output specie 1 
	output= open(path+specie1+'_BBHrG4.txt',"w") # file opening for reading
	output.write("\n".join(BBHG4sSpecie1))
	output.close()
						
	# write output specie 1 
	output= open(path+specie2+'_BBHrG4.txt',"w") # file opening for reading
	output.write("\n".join(BBHG4sSpecie2))
	output.close()
#----------------------------------------------------------------------#	
main()


	










