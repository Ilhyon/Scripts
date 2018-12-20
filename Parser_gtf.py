#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
from pprint import pprint


#----------------------------------------------------------------------#
def ImportFeature(filename):
	""" 
		Parameters
	    ----------
	     filename : string, name of the file containing feature to parse
	    Returns
	    -------
	"""
	dicoFeature = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: # browse all lines
			words=l.split('\t')
			attribute = words[9].split(';')
			idGene = attribute[1].split('"')[2];
			if idGene not in dicoFeature :
				chrm = word[1]
				feature = words[3]
				startFeature = words[4]
				endFeature = words[5]
				strand = words[7]
				biotype = 
	return(DicoBBH)
