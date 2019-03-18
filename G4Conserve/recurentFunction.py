#!/usr/bin/env python
# -*- coding: utf-8 -*-:

def createDictionaryBiotypeByTranscript(filename):
	""" Create dictionary with as value the biotype of the transcript 
		for all transcripts of the chromosome

	    Parameters
	    ----------
	    filename : string
		name of file which contain the liste of biotype for each 
		transcript of this chromosome

	    Returns
	    -------
	    dico
		dictionary {IdTr : biotype}
	"""
	dico = {}
	inputfile = open(filename,'r')
	for line in inputfile:
		words = line.split('|')
		idTr = words[1]
		biotypeTranscript = words[3].rstrip()
		if idTr not in dico:
			dico[idTr] = biotypeTranscript
	inputfile.close()
	return dico

def GetAnnotationTranscript(filename, Coding, BiotypeByTranscript):
	""" Return a dictionnary which define if the annotation of a
		transcript is good.
	    Parameters
	    ----------
	    filename : string
		name of the file which contain informatiob fr each transcript
		(except the transcript's biotype)
	    Coding : ensemble of string, in list
		name of biotype which form the suprafamily protein coding
	    BiotypeByTranscript : dico
		dictionary, where key =id of transcript and value = biotype of this transcript
	    Returns
	    -------
	    AnnotationTranscript: dico {idTr : bool}
	"""
	AnnotationTranscript = {}
	inputfile = open(filename,"r")
	for line in inputfile:
		words = line.split('|')
		transcriptId = words[0]
		start5 = words[7]
		end5 = words[8]
		start3 = words[9] 
		end3 = words[10].rstrip() 
		answer = True # variable answer by defaul True
		transcriptBiotype = ''
		transcriptBiotype = BiotypeByTranscript.get(transcriptId)
		if transcriptBiotype not in Coding:
			if (start5!='' or end5!='' or start3!='' or end3!='' ):
				# but if transcript has a 5'UTR or an 3' UTR bad anotation
				answer = False
		AnnotationTranscript[transcriptId] = answer
	inputfile.close()
	return AnnotationTranscript
