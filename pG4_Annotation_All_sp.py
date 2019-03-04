#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import sys
import os
import scipy.sparse as sparse
import re
import argparse
from pprint import pprint

def mean(liste): # compute the mean of a list
   	return sum(liste)/len(liste)

def updateEndG4(strand, windowEnd, windowStart):
	"""
		when a new window is add to a pG4
	"""
	if strand == "1":
		endG4 = windowEnd
	else:
		endG4 = windowStart
	return endG4

def retrieveG4CoordByStrand(strand,
						  startWindow,
						  endWindow):
	if strand == "1":
		startG4 = startWindow
		endG4 = endWindow
	else :
		startG4 = endWindow
		endG4 = startWindow
	return startG4, endG4

def createIdG4(gene, startG4, endG4, strand):
	if strand == "1":
		idG4 = gene+"|"+str(startG4)+"|"+str(endG4)+"|"+strand
	elif strand == "-1":
		idG4 = gene+"|"+str(endG4)+"|"+str(startG4)+"|"+strand
	return idG4

def addWindowToSequenceG4(g4Seq, windowSeq, step, windowLength):
	if len(windowSeq) < windowLength: 
		# if the sequence is shorter than the sequence of the windows
		# (generraly last windows above the thresolds)
		g4Seq += windowSeq[-(len(windowSeq)-(windowLength-step)):]
	else: #
		g4Seq += windowSeq[-step:] # take the stepsize added
	return g4Seq

def positionChromosomiqueGenePositif(position,
									junLength,
									startIntron,
									endIntron): # convert a position of an Intron in chromosomique position (gene +)
	""" 
		Give chromosomal positions of G4 around a junction,
		in a forward gene.
	"""
	if position <= junLength:	# because start G4 classifier from 0
		position = startIntron - 1 - junLength + position
	else:	# downstream sequence
		position = endIntron + 1 - junLength + position - 1
	return position

def positionChromosomiqueGeneNegatif(position,
									junLength,
									startIntron,
									endIntron):	# convert a position of an Intron in chromosomique position (gene -)
	""" 
		Give chromosomal positions of G4 around a junction,
		in a reverse gene.
	"""
	if position <= junLength + 1:	# because strart G4 classifier from 0
		position = startIntron + 1 + junLength + 1 - position
	else:	# if from sequence aval
		position = endIntron + junLength - position
	return position

def chromosomalPosition(position, 
						extension,
						startIntron,
						endIntron,
						strand): # convert a position of an Intron in chromosomique position (gene +)
	""" 
		Give real position (chromosomique) of start or end of G4 around 
		a junction from coordination of an intron FOR A GENE POSITIF
	"""
	if strand == "1" :
		if (int(position) <= extension):	## because start G4 classifier from 0
			position=int(startIntron)-1-extension+int(position)
		else:	# if from downstream sequence
			position=int(endIntron)+1-extension+int(position)-1
		return position
	else :
		if (int(position) <= extension+1):	## because strart G4 classifier from 0
			position=int(startIntron)+1+extension+1-int(position)
		else:	# if from downstream sequence
			position=int(endIntron)+extension-int(position)
		return position

def CreateDictionaryBiotypeByTranscript (filename):
	""" 
		Create dictionary with as value the biotype of the transcript 
		for all transcripts of the chromosome
	"""
	dico={}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for line in lines: # for each transcript
			if line :
				words=line.split('|')
				transcript=words[1]
				biotypeTranscript=words[3]
				if transcript not in dico :
					dico[transcript]=biotypeTranscript
	return dico

def CreateDictionaryStrandByGene (filename): 
	""" 
		Create dictionary with as value the strand for all genes 
		of the chromosome.
	"""
	dico={}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for line in lines: # for each transcript
			if line :
				words=line.split('|')
				gene=words[1].rstrip()
				strand=words[3].rstrip()
				if gene not in dico:
					dico[gene]=strand
	return dico

def readLineG4Screener(line, StrandByGene, feature):
	"""
		There is many information in the output of G4RNA screener so 
		they will be contain in a dictionary for more readability.
		The last argument is the feature type (gene or junction).
		This allow to get "good" names for keys in dico.
	"""
	words = line.split('\t') # get all element in a list
	dicoLine = {feature+"Description" : words[1],
				"GeneID" : words[1].split("|")[0],
				"Strand" : StrandByGene[words[1].split("|")[0]],
				feature+"Start" : int(words[1].split("|")[1]),
				feature+"End" : int(words[1].split("|")[2]),
				"cGcC" : float(words[2]),
				"g4H" : float(words[3]),
				"WindowSeq" : words[4],
				"WindowStart" : int(words[5]),
				"WindowEnd" : int(words[6]),
				"g4NN" : float(words[7])}
	return dicoLine

def windowOverUpdateDicoTmppG4(dicoTmp, dicoLine,
								startG4, endG4,
								dicoParametersTool,
								feature):
	"""
		Storage of the information for window over all the threshold.
	"""
	if startG4 not in dicoTmp:
		dicoTmp.update({startG4 : {"cGcC" : [dicoLine["cGcC"]],
									"g4H" : [dicoLine["g4H"]],
									"g4NN" : [dicoLine["g4NN"]],
									"End" : endG4,
									"g4Sequence" : dicoLine["WindowSeq"],
									feature+"Description" : dicoLine[feature+"Description"]}})
	else:
		dicoTmp[startG4]["cGcC"].append(dicoLine["cGcC"])
		dicoTmp[startG4]["g4H"].append(dicoLine["g4H"])
		dicoTmp[startG4]["g4NN"].append(dicoLine["g4NN"])
		dicoTmp[startG4]["End"] = endG4
		g4Seq = addWindowToSequenceG4(dicoTmp[startG4]["g4Sequence"],
									dicoLine["WindowSeq"],
									dicoParametersTool["step"],
									dicoParametersTool["windowLength"])
		dicoTmp[startG4]["g4Sequence"] = g4Seq
	return dicoTmp

def windowUnderUpdateDicoTmppG4(dicoTmp, dicoLine, G4DetectedInGene,
								startG4):
	"""
		Processing of the information if we are under a threshold
	"""
	if startG4 in dicoTmp :
		meanCGcC = mean(dicoTmp[startG4]["cGcC"])
		meanG4Hunter = mean(dicoTmp[startG4]["g4H"])
		meanG4NN = mean(dicoTmp[startG4]["g4NN"])
		idG4 = createIdG4(dicoLine["GeneID"],
						startG4, dicoTmp[startG4]["End"],
						dicoLine["Strand"])
		if idG4 not in G4DetectedInGene and dicoLine["Strand"]: 
			G4DetectedInGene[idG4] = {"cGcC" : str(meanCGcC),
									"g4H" : str(meanG4Hunter),
									"pG4Sequence" : dicoTmp[startG4]["g4Sequence"],
									"g4NN" : str(meanG4NN)}
			del dicoTmp[startG4]			
	return G4DetectedInGene, dicoTmp

def ReturnG4InGene(G4DetectedInGene,
					inputfile,
					dicoParametersTool,
					StrandByGene,
					feature):
	"""
		Return a dictionary with informations of region G4 detected 
		(score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		regions of G4 discovered in genes.
	"""
	dicoTmp = {}
	windowCounter = 0
	startG4 = None
	descriptionGeneOverThreshold = ""
	with open(inputfile) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for line in lines:
			if (re.search('^[0-9]', line)): # if line is not the header
				dicoLine = readLineG4Screener(line, StrandByGene, feature)
				if descriptionGeneOverThreshold:
					if descriptionGeneOverThreshold != dicoLine[feature+"Description"]:
						# before checking scores, we need to know if we are
						# still in the same gene
						windowCounter = 0
				if (dicoLine["cGcC"] >= dicoParametersTool["cGcC"] and
					dicoLine["g4H"] >= dicoParametersTool["g4H"] and
					dicoLine["g4NN"] >= dicoParametersTool["g4NN"]):
					# all scores are over the threshold
					windowCounter += 1
					descriptionGeneOverThreshold = dicoLine[feature+"Description"]
					if windowCounter == 1:
						startG4, endG4 = retrieveG4CoordByStrand(dicoLine["Strand"],
											dicoLine["WindowStart"],
											dicoLine["WindowEnd"])
					elif windowCounter > 1:
						endG4 = updateEndG4(dicoLine["Strand"],
									dicoLine["WindowEnd"],
									dicoLine["WindowStart"])
					dicoTmp = windowOverUpdateDicoTmppG4(dicoTmp,
														dicoLine,
														startG4,
														endG4,
														dicoParametersTool,
														feature)
				elif (dicoLine["cGcC"] < dicoParametersTool["cGcC"] or
					dicoLine["g4H"] < dicoParametersTool["g4H"] or
					dicoLine["g4NN"] < dicoParametersTool["g4NN"]):
					# one of the score is under the threshold or 
					# this windows contain info of an other gene
					if startG4:
						G4DetectedInGene, dicoTmp = windowUnderUpdateDicoTmppG4(dicoTmp,
																			dicoLine,
																			G4DetectedInGene,
																			startG4)
						windowCounter = 0
					descriptionGeneOverThreshold = dicoLine[feature+"Description"]
	return G4DetectedInGene

def G4IsOnJunction(startG4, endG4, startBorder, endBorder):
	""" 
		Return a boolean : if the G4 is on the junction it's true 
		elswise it's false.
	"""
	onJunction=False
	if (startG4 <startBorder and endG4 >startBorder) or (startG4 >startBorder and endG4<startBorder):
		onJunction=True
	return onJunction

def ReturnG4InJunction(inputfile,
					dicoParam,
					StrandByGene,
					feature): ## for G4 in junction 
	"""
		Return a dictionary with informations of region G4 detected 
		(score cgCc, scote G4Hunter, sequence and score G4NN) for each 
		regions of G4 discovered in junctions.
	"""
	dicoTmp = {}
	windowCounter = 0
	startG4 = None
	descriptionGeneOverThreshold = ""
	with open(inputfile) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for line in lines:
			if (re.search('^[0-9]', line)): # if the line is not the header of the file 
				dicoLine = readLineG4Screener(line, StrandByGene,"Ah")
				if (cGcC >= dicoParam["cGcC"] and
					g4H >= dicoParam["g4H"] and
					g4NN >= dicoParam["g4NN"] and
					startWindow):
					# if window contain value > several threshold
					onJunction = False
					passed = True # passage over the differents threshold for this window
					if oldPassed != passed: # if it's the first windows (beginning if the passage)
						listeCGcC = [] # list which will contain the values of score cgCG for windows over the THRESHOLD_CGCC
						listeG4H = []  # list which will contain the values of score G4H for windows over the THRESHOLD_G4H
						listeG4NN = []  # list which will contain the values of score G4NN for windows over the THRESHOLD_G4NN
						descriptionOverThreshold = description # assignation of description for this variable
						gene_startG4 = gene
						sequenceG4 = sequenceWindow
						startFirstWindow = startBorder
						endFirstWindow = endBorder
						oldPassed = passed # assignation, we pass over the thresholds
						#if (strand == str(1)): # if gene positif
						startG4 = int(startWindow) # startG4 same as startWindow of G4 screener
						endG4 = int(endWindow) # endG4 same as endWindow of G4 screener
						#else:	# if gene negatif
							#startG4=int(endBorder)-(int(startWindow)-int(startBorder)) # calcul of startG4
							#endG4=int(startG4)-(len(sequence))+1 # calcul of endG4
						listeCGcC = [cGcC] # add score cGCc  in the list 
						listeG4Hunter = [g4H] # add score G4H  in the list 
						listeG4NN = [g4NN] # add score G4NN  in the list 	
					else: # if it's not the first windows above the thresholds
						if len(sequenceWindow) < dicoParam["windowLength"]: 
						# if the sequence is shorter than the window's 
						# length (can happen for the last one)
							overlap = dicoParam["windowLength"]
							- dicoParam["step"]
							sequenceG4 += sequenceWindow[-len(sequenceWindow) - overlap :]
							# we only retrieve the end of the sequence
						else:
							sequenceG4 += sequenceWindow[- dicoParam["step"]:] # take the stepsize added
							# we only retrieve the sequence that don't
							# overlap the previous window
						#if (strand == str(1)): # if gene positif
						endG4=endWindow # endG4 moves from the end of the first window,  same as endWindow of G4 screener
						#else: # if gene negatif
							#endG4=int(startG4)-(len(sequenceG4))+1 # endG4 moves from the end of the first window, # calcul of endG4
						lastRow=numRow # numRow of the last windows study over the threshols
						listeCGcC.append(cGcC) # add score cGcC for the windows over the threshold
						listeG4Hunter.append(g4H) # add score G4H for the windows over the threshold
						listeG4NN.append(g4NN) # add score G4NN for the windows over the threshold
					
				if (cGcC < dicoParam["cGcC"] or
					g4H < dicoParam["g4H"] or
					g4NN < dicoParam["g4NN"] or
					descriptionOverThreshold != description):
				# if one of the score is under his threshold or if this windows contain info of an other gene
					passed=False # passed passed became false (pass under the thresolds)
					if (oldPassed != passed ): # last windows before under the thresolds
						meanCGcC=mean(listeCGcC) # mean of the score cGcC will be the score for this region containing G4
						meanG4Hunter=mean(listeG4Hunter) # mean of the score G4Hunter will be the score for this region containing G4
						meanG4NN=mean(listeG4NN) # mean of the score G4NN will be the score for this region containing G4
						oldPassed=passed # assignation, we pass under the thresholds
						if (strand == "1"): # if gene positif
							startG4=positionChromosomiqueGenePositif(startG4, dicoParam["junctionLength"], startFirstWindow, endFirstWindow)	
							endG4=positionChromosomiqueGenePositif(endG4, dicoParam["junctionLength"], startFirstWindow, endFirstWindow)
						else: # if gene negatif
							startG4=positionChromosomiqueGeneNegatif(startG4, dicoParam["junctionLength"], startFirstWindow, endFirstWindow)	
							endG4=positionChromosomiqueGeneNegatif(endG4, dicoParam["junctionLength"], startFirstWindow, endFirstWindow)
						if (strand == "1"):
							headerG4=gene+"|"+str(startG4)+"|"+str(endG4)+"|"+strand
						elif (strand == "-1"):
							headerG4=gene+"|"+str(endG4)+"|"+str(startG4)+"|"+strand
						onJonction=G4IsOnJunction(startG4, endG4, startBorder, endBorder)	
						if headerG4 not in G4DetectedInJunction and  onJonction :
							G4DetectedInJunction[headerG4]=str(meanCGcC), str(meanG4Hunter),sequenceG4 , str(meanG4NN)
	return G4DetectedInJunction	

def BorderOfTranscript(start5 , end5 , start3 , end3, exonList, intronList, strand):
	listBorder=[]
	if not start5 :	# if transcript doesn't have 5'
		if not end3: # if transcript doesn't have 3'
			# border will be given by exons and introns
			for couple_exon in exonList:	# for each couple of exon of this transcript
				if couple_exon :
					startExon=couple_exon.split("-")[0] # start of this exon
					endExon=couple_exon.split("-")[1] # start of this exon
					listBorder.append(int(startExon))	# add information in list
					listBorder.append(int(endExon))	# add information in list
			for couple_Intron in intronList: # for each couple of intron of this transcript
				if couple_Intron :
					startIntron=couple_Intron.split("-")[0]	# start of this intron
					endIntron=couple_Intron.split("-")[1]	# end of this intron
					listBorder.append(startIntron)	# add information in list
					listBorder.append(endIntron)	# add information in list
			if strand == "1":	#if gene positif, 5'< 3'
				borderInferior=min(listBorder)
				borderSuperior=max(listBorder)
			else: 	#if gene negatif, 5'> 3'
				borderInferior=max(listBorder)
				borderSuperior=min(listBorder)
		else: # if transcript doesn't have 5'
			# border will be given by exons and introns, and 3'
			for couple_exon in exonList:	# for each couple of exon of this transcript
				if couple_exon :
					startExon=couple_exon.split("-")[0] # start of this exon
					endExon=couple_exon.split("-")[1] # start of this exon
					listBorder.append(int(startExon))	# add information in list
					listBorder.append(int(endExon))	# add information in list
			for couple_Intron in intronList: # for each couple of intron of this transcript
				if couple_Intron :
					startIntron=couple_Intron.split("-")[0]	# start of this intron
					endIntron=couple_Intron.split("-")[1]	# end of this intron
					listBorder.append(int(startIntron))	# add information in list
					listBorder.append(int(endIntron))	# add information in list
			if (strand == "1"): # if gene positif
				borderInferior=min(listBorder)
			else:	# if gene negatif
				borderInferior=max(listBorder)
			borderSuperior=int(end3)
	elif not start3 : ##  if transcript doesn't have 3'
		# border will be given by exons and introns and 5'
		listBorder.append(start5)
		listBorder.append(end5)
		for couple_exon in exonList:
			if couple_exon :
				startExon=couple_exon.split("-")[0]
				endExon=couple_exon.split("-")[1]
				listBorder.append(int(startExon))
				listBorder.append(int(endExon))
		for couple_Intron in intronList:
			if not (couple_Intron==""):
				startIntron=couple_Intron.split("-")[0]
				endIntron=couple_Intron.split("-")[1]
				listBorder.append(int(startIntron))
				listBorder.append(int(endIntron))
		if (strand == "1"):
			borderSuperior=max(listBorder)
		else:
			borderSuperior=min(listBorder)
		borderInferior=int(start5)
	else: ##  if transcript has 5', 3', exons and introns
		borderInferior=int(start5)
		borderSuperior=int(end3)
	return [borderInferior,borderSuperior]

def G4IsInTranscript(strand, coordG4, borderTranscript):	
	""" 
		Give boolean, if a G4 is located in a transcript
	"""
	startG4=int(coordG4[0])
	endG4=int(coordG4[1])
	borderInferior=int(borderTranscript[0])
	borderSuperior=int(borderTranscript[1])
	inTranscript=False
	if strand == "1" :
		if (int(startG4) >= int(borderInferior) and int(startG4) <= int(borderSuperior) and int(endG4) >= int(borderInferior) and int(endG4) <= int(borderSuperior)):
			inTranscript=True
	else :
		if (int(startG4) <= int(borderInferior) and int(startG4) >= int(borderSuperior) and int(endG4) <= int(borderInferior) and int(endG4) >= int(borderSuperior)):
			inTranscript=True
	return inTranscript

def GetLocalisationPositif(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList):
	""" 
		Give localisation of a G4 region in a forward transcript
	"""
	startG4=int(coordG4[0]) 
	endG4=int(coordG4[1])
	borderInferior=int(borderTranscript[0])
	borderSuperior=int(borderTranscript[1]) 
	localisation='NAP'
	if transcriptId in BiotypeByTranscript :
		biotypeTranscript = BiotypeByTranscript[transcriptId]
	else:
		biotypeTranscript = "No biotype"
	if biotypeTranscript in ProteinCoding:
		codant='CDS'
	else:
		codant='ExonNC'
	if end5 and startG4 <= int(end5) :	# if G4 is positioned in the 5' region
		if (endG4 <= int(end5)): # if G4 is positioned totally in the 5' region
			localisation="5"
		else:	# if G4 is positioned in the 5' region with overlapping
			localisation="junction_5_"+codant
	elif start3 and endG4 >= int(start3) :	# if G4 is positioned in the 3' region
		if (startG4 >= int(start3)):	 # if G4 is positioned totally in the 3' region
			localisation="3"
		else:	# if G4 is positioned in the 3' region with overlapping
			localisation="junction_"+codant+"_3" 
	else:
		for couple_exon in exonList: # if G4 is positioned in the exonic region, for each exon
			if couple_exon : # if transcript contain exon
				startExon=int(couple_exon.split("-")[0])
				endExon=int(couple_exon.split("-")[1])
				if (startG4 >= startExon and startG4 <= endExon):
					if (endG4 <= endExon):
						localisation=codant
					else:
						if (codant =='CDS'):
							localisation="junction_"+codant+"_Intron"
						else:
							localisation="junction_"+codant+"_IntronNC"
		for couple_Intron in intronList:	# if G4 is positioned in the intronic region, for each intron
			if couple_Intron :	# if transcript contain intron
				startIntron=int(couple_Intron.split("-")[0])
				endIntron=int(couple_Intron.split("-")[1])
				if (startG4 >= startIntron and startG4 <= endIntron):
					if (endG4 <= endIntron):
						if (codant=='CDS'):
							localisation="Intron"
						else:
							localisation="IntronNC"
					else:
						if (codant=='CDS'):
							localisation="junction_Intron_"+codant
						else:
							localisation="junction_IntronNC_"+codant
	if localisation or localisation != 'NA' :
		return localisation

def GetLocalisationNegatif(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList ):
	""" 
		Give the localisation of a G4 region in the reverse transcript
	"""
	startG4=int(coordG4[0])
	endG4=int(coordG4[1])
	borderInferior=int(borderTranscript[0])
	borderSuperior=int(borderTranscript[1]) 
	localisation='NAN'
	if transcriptId in BiotypeByTranscript :
		biotypeTranscript = BiotypeByTranscript[transcriptId]
	else:
		biotypeTranscript = "No biotype"
	if (biotypeTranscript in ProteinCoding):
		codant='CDS'
	else:
		codant='ExonNC'
	if end5 and endG4 >= int(end5):	# if G4 is positioned in the 5' region
		if (startG4 >= int(end5)): # if G4 is positioned totally in the 5' region
			localisation="5"
		else:	# if G4 is positioned in the 5' region with overlapping
			localisation="junction_5_"+codant
	elif start3 and startG4 <= int(start3) :	# if G4 is positioned in the 3' region
		if (endG4 <= int(start3)):	 # if G4 is positioned totally in the 3' region
			localisation="3"
		else:	# if G4 is positioned in the 3' region with overlapping
			localisation="junction_"+codant+"_3" 
	else:
		for couple_exon in exonList: # if G4 is positioned in the exonic region, for each exon
			if couple_exon : # if transcript contain exon
				startExon=int(couple_exon.split("-")[0])
				endExon=int(couple_exon.split("-")[1])
				if (endG4 <= startExon and endG4 >= endExon):
					if (startG4 >= endExon):
						localisation=codant
					else:
						if (codant =='CDS'):
							localisation="junction_"+codant+"_Intron"
						else:
							localisation="junction_"+codant+"_IntronNC"
		for couple_Intron in intronList:	# if G4 is positioned in the intronic region, for each intron
			if couple_Intron :	# if transcript contain intron
				startIntron=int(couple_Intron.split("-")[0])
				endIntron=int(couple_Intron.split("-")[1])
				if (endG4 <= startIntron and endG4 >= endIntron):
					if (startG4 >= endIntron):
						if (codant=='CDS'):
							localisation="Intron"
						else:
							localisation="IntronNC"
					else:
						if (codant=='CDS'):
							localisation="junction_Intron_"+codant
						else:
							localisation="junction_IntronNC_"+codant
	if localisation  or localisation != 'NA' :
		return localisation

def GetAnnotationTranscript(filename,ProteinCoding,BiotypeByTranscript):
	""" 
		This fonction define if the annotation of a transcript is good.
	"""
	AnnotationTranscript={} # creation dico empty
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for line in lines:
			if line :
				words = line.split('|')
				transcriptId = words[0].rstrip() 
				start5 = words[7].rstrip()
				end5 = words[8].rstrip() 
				start3 = words[9].rstrip() 
				end3 = words[10].rstrip() 
				answer = True # variable answer by defaul True
				if transcriptId in BiotypeByTranscript :
					biotypeTranscript = BiotypeByTranscript[transcriptId]	
					# get transcriptBiotype of this TRanscriptId
					if biotypeTranscript not in ProteinCoding: 
						# if transcript not a protein coding
						if start5 != '' or end5 != '' or start3 != '' or end3 != '' : 
							# but if transcript has a 5'UTR or an 3' UTR
							answer = False # has not a good annotation in Ensembl
					AnnotationTranscript[transcriptId] = answer
				else:
					AnnotationTranscript[transcriptId] = None		
	return AnnotationTranscript

def GetLocalisationG4InJunction(BiotypeByTranscript,ProteinCoding,transcriptId):
	localisation='NAJ'
	biotypeTranscript= BiotypeByTranscript[transcriptId]
	if biotypeTranscript in ProteinCoding :
		codant='CDS'
	else:
		codant='ExonNC'
	return "junction_"+codant+"_"+codant

def GetLocalisationG4InTranscript(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList, strand):
	""" 
		Define which fonction will be used for determinate the 
		localisation in transcript. 
	"""
	if strand ==  "1" :
		localisation=GetLocalisationPositif(BiotypeByTranscript,
											ProteinCoding,
											transcriptId,
											borderTranscript,
											coordG4,
											end5,
											start3,
											exonList,
											intronList)
	else :
		localisation=GetLocalisationNegatif(BiotypeByTranscript,
											ProteinCoding,
											transcriptId,
											borderTranscript,
											coordG4,
											end5,
											start3,
											exonList,
											intronList)
	return localisation

def AddG4InTranscriptome(G4InTranscript,transcriptId, descriptionG4,informationsOfG4,localisationInTranscript, biotypeTranscript):
	""" 
		Add for each G4 predicted : the informations from G4 screener 
		(cGcC, g4H, sequence and G4NN), the localisation of G4 and the 
		biotype of transcript
	"""
	value = [] # list empty of value for a g4 by gene
	for information in informationsOfG4: # add information already obtained by G4 screener (cGcC, g4H, sequence and G4NN)
		value.append(information)
	value.append(localisationInTranscript) # add localisation of this G4
	value.append(biotypeTranscript) # add biotype of this G4 for this transcript
	headerG4=transcriptId+'|'+descriptionG4 # header uniq by G4 in transcript
	if headerG4 not in G4InTranscript : # if the transcript doesn't be traited before
		G4InTranscript[headerG4]=value # add this G4 in the dictionary
	return G4InTranscript	# return the dictionary

def AddG4InGenome(G4InGenome, geneId, descriptionG4, localisationInTranscript):
	""" 
		Add for each G4 predicted in the gene all localisation possible
	"""
	headerG4=geneId+':'+descriptionG4  # header uniq by G4 in gene
	if headerG4 not in G4InGenome : # if the gene doesn't be traited before
		G4InGenome[headerG4]=localisationInTranscript 	# traited it
	else:	# if the gene already traited before
		localisationInGene=G4InGenome[headerG4] # get the localisation already known for this G4
		if localisationInGene and localisationInTranscript not in localisationInGene:
			G4InGenome[headerG4]=localisationInGene+';'+localisationInTranscript # add this localisation
	return G4InGenome

def AddTranscriptPerG4(TranscriptPerG4, descriptionG4,transcriptId): 
	if descriptionG4 not in TranscriptPerG4 : 
		TranscriptPerG4[descriptionG4] = transcriptId
	else:
		if TranscriptPerG4 : # id dico not empty
			lastTranscript= TranscriptPerG4[descriptionG4]
			TranscriptPerG4[descriptionG4] = lastTranscript +"-"+transcriptId
	return TranscriptPerG4

def JuncionIsInTranscript(coordG4,intronList):
	startG4=int(coordG4[0]) 
	endG4=int(coordG4[1])
	inTranscript=False
	for couple_Intron in intronList: # for each couple of intron of this transcript
		if couple_Intron :
			endIntron=int(couple_Intron.split("-")[1])	# end of this intron
			if startG4 <= endIntron and endIntron<= endG4 :
				inTranscript=True
	return inTranscript

def GetlisteG4InGene(G4Detected, listeG4InGene):
	""" 
		Create a dictionary with each gene and the list of G4 predicted
		in it.
	"""
	for descriptionG4 in G4Detected:
		geneID=descriptionG4.split('|')[0]
		if geneID not in listeG4InGene :	
			listeG4=[]
		else:
			listeG4=listeG4InGene[geneID]
		listeG4.append(descriptionG4)
		listeG4InGene[geneID]=listeG4
	return listeG4InGene

def ExtractionG4InTranscript(directory, specie, chromosome, G4InTranscript):
	output= open(directory+"/"+specie+"_chr"+chromosome+"_G4InTranscript.txt","w") ## file opening
	output.write("InfoG4ByTranscript\tcGcC\tG4Hunter\tsequenceG4\tG4NN\tlocalisation\ttranscriptBiotype\n")
	for G4 in G4InTranscript :
		if None not in G4InTranscript[G4]: # because some transcriptID from ensembl donMt contain info of biotype (as ENST00000604369)
			output.write(G4+"\t"+'\t'.join(G4InTranscript[G4])+"\n")
	output.close()

def ExtractionG4InGenome(directory, specie, chromosome, G4InGenome):
	output= open(directory+"/"+specie+"_chr"+chromosome+"_G4InGenome.txt","w") ## file opening
	output.write("InfoG4ByGene\tLocalisation(s)\n")
	for G4 in G4InGenome :
		output.write(G4+"\t"+G4InGenome[G4]+"\n")
	output.close()

def ExtractionTranscriptPerG4(directory, specie, chromosome, TranscriptPerG4):
	output= open(directory+"/"+specie+"_chr"+chromosome+"_TranscriptPerG4.txt","w") ## file opening
	output.write("CoordonnéesG4\tTranscript(s)\n")
	for G4 in TranscriptPerG4 :
		output.write(G4+"\t"+TranscriptPerG4[G4]+"\n")
	output.close()

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation : filter the output of G4 screener and add some informations')
	GITDIR=os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR+'data')
	parser.add_argument ('-chr', '--chromosome', default = 'X')
	parser.add_argument ('-sp', '--specie', default = 'HS')
	parser.add_argument ('-g4h', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-cgcc', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-g4nn', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-e', '--EXTENSION', default = 100)
	parser.add_argument ('-w', '--WINDOW', default = 60)
	parser.add_argument ('-s', '--STEP', default = 10)
	return parser

def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path # directory which contain all the directory chromosome
	chromosome=arg.chromosome
	specie=arg.specie
	dicoParameterG4Screener = {"g4H" : float(arg.THRESHOLD_G4H),
							"cGcC" : float(arg.THRESHOLD_CGCC),
							"g4NN" : float(arg.THRESHOLD_G4NN),
							"junctionLength" : int(arg.EXTENSION),
							"windowLength" : int(arg.EXTENSION),
							"step" : int(arg.STEP)}

	GITDIR=os.getcwd()
	# ~ print "Chromosome : "+chromosome
	ProteinCoding=['IG_C_gene',
					'IG_D_gene',
					'IG_J_gene',
					'IG_LV_gene',
					'IG_M_gene',
					'IG_V_gene',
					'IG_Z_gene',
					'nonsense_mediated_decay',
					'nontranslating_CDS',
					'non_stop_decay',
					'protein_coding',
					'TR_C_gene',
					'TR_D_gene',
					'TR_gene',
					'TR_J_gene',
					'TR_V_gene']
	
	G4DetectedInGene = {}
	listeG4InGeneEntire={}
	listeG4InGeneJunction={}

	directory=path+'/chr'+chromosome # variable directory which contain the data for this chromosome
	index=directory+'/'+specie+'_transcript_unspliced_chr'+chromosome+'_Index.txt'	# file which contain info by transcript for this chromosome
	indexBiotypeTranscript=path+'/transcriptType/transcriptType_chr'+chromosome	# file which contain biotype of transcript for this chromosome

	BiotypeByTranscript = CreateDictionaryBiotypeByTranscript(indexBiotypeTranscript)
	StrandByGene = CreateDictionaryStrandByGene(index)
	
	AnnotationTranscript=GetAnnotationTranscript(index,ProteinCoding,BiotypeByTranscript) # annotation transcript
		
	#####Creation dictionary G4_detected_gene and G4_detected_junction , for G4 found in gene and junction with scores > thresolds
	for path, dirs, files in os.walk(directory): # for each element of the directory to passed
		for filename in files: # for each files
			inputfile=directory+'/'+filename
			##EXTRACTION OF G4 BY SCORE G4NN > THRESHOLD IN A DICO
			if ('gene_unspliced' in filename and '.csv' in filename ): ## for G4 in gene	
				G4DetectedInGene=ReturnG4InGene(G4DetectedInGene,
												inputfile,
												dicoParameterG4Screener,
												StrandByGene,
												"Gene")
				print len(G4DetectedInGene)
			elif ('transcript_unspliced' in filename and '.csv' in filename): ## for G4 in junction CDS-CDS --> from splicing
				G4DetectedInJunction=ReturnG4InJunction(inputfile,
										dicoParameterG4Screener,
										StrandByGene)
	
	listeG4InGeneEntire=GetlisteG4InGene(G4DetectedInGene, listeG4InGeneEntire)
	listeG4InGeneJunction=GetlisteG4InGene(G4DetectedInJunction, listeG4InGeneJunction)
	
	G4InTranscript={}
	G4InGenome={}
	TranscriptPerG4={}
	
############################################################################## g4 entire
	inputfile= open(index,"r") # file opening for reading
	for line in inputfile: # for each transcript
		words=line.split('|')
		transcriptId=words[0].rstrip() 
		geneId=words[1].rstrip() 
		chromosome=words[2].rstrip() 
		strand=words[3].rstrip() 
		biotypeGene=words[4].rstrip() 
		exonList=words[5].rstrip().split(";") ## transform in array
		intronList=words[6].rstrip().split(";") ## transform in array
		### resolution probleme start5 and end5 for gene negatif
		if (strand == str(1)):
			start5=words[7].rstrip()
			end5=words[8].rstrip()
			start3=words[9].rstrip()
			end3=words[10].rstrip()
		else:
			end5=words[7].rstrip()
			start5=words[8].rstrip() 
			end3=words[9].rstrip()
			start3=words[10].rstrip()
		if transcriptId in BiotypeByTranscript:
			biotypeTranscript=BiotypeByTranscript[transcriptId] # get biotype of this transcript
		else:
			biotypeTranscript = "No biotype"
		borderTranscript=BorderOfTranscript(start5,
											end5,
											start3,
											end3,
											exonList,
											intronList,
											strand)
		if geneId in listeG4InGeneEntire : # if no G4 predicted in that gene
			listeG4InGene = listeG4InGeneEntire[geneId]
			for G4InGene in listeG4InGene:
				startG4 = int(G4InGene.split('|')[1])
				endG4 = int(G4InGene.split('|')[2])
				coordG4 = [startG4,endG4]
				g4inTranscript = G4IsInTranscript(strand, coordG4, borderTranscript)
				annotationTranscript = AnnotationTranscript[transcriptId]
				if g4inTranscript and annotationTranscript :
					informationsOfG4 = list(G4DetectedInGene[G4InGene])
					localisationInTranscript = GetLocalisationG4InTranscript(BiotypeByTranscript,ProteinCoding,transcriptId,borderTranscript,coordG4, end5, start3, exonList, intronList, strand)	
					descriptionG4 = chromosome+':'+str(startG4)+'-'+str(endG4)+'|'+strand
					G4InTranscript = AddG4InTranscriptome(G4InTranscript,transcriptId, descriptionG4,informationsOfG4,localisationInTranscript, biotypeTranscript)
					G4InGenome = AddG4InGenome(G4InGenome, geneId, descriptionG4, localisationInTranscript)
					TranscriptPerG4 = AddTranscriptPerG4(TranscriptPerG4, descriptionG4,transcriptId)

	############################################################################## g4 junction
	inputfile = open(index,"r") # file opening for reading
	for line in inputfile: # for each transcript
		words = line.split('|')
		transcriptId = words[0] 
		geneId = words[1] 
		chromosome = words[2] 
		strand = words[3] 
		biotypeGene = words[4] 
		exonList = words[5].split(";") ## transform in array
		intronList = words[6].split(";") ## transform in array
		### resolution probleme start5 and end5 for gene negatif
		if strand == "1":
			start5 = words[7]
			end5 = words[8]
			start3 = words[9]
			end3 = words[10]
		else:
			end5 = words[7]
			start5 = words[8] 
			end3 = words[9]
			start3 = words[10]
		if transcriptId in BiotypeByTranscript:
			biotypeTranscript = BiotypeByTranscript[transcriptId]
		else:
			biotypeTranscript = "No biotype"
		if geneId in listeG4InGeneJunction:
			listeG4InGene=listeG4InGeneJunction[geneId]
			for G4InGene in listeG4InGene:
				startG4=int(G4InGene.split("|")[1])
				endG4=int(G4InGene.split("|")[2])
				coordG4=[startG4,endG4]
				junctionInTranscript=JuncionIsInTranscript(coordG4,intronList)
				annotationTranscript=AnnotationTranscript[transcriptId]
				if junctionInTranscript and annotationTranscript :
					informationsOfG4=list(G4DetectedInJunction[G4InGene])
					localisationInTranscript=GetLocalisationG4InJunction(BiotypeByTranscript,ProteinCoding,transcriptId)
				
					descriptionG4=chromosome+':'+str(startG4)+'-'+str(endG4)+'|'+strand
					G4InTranscript=AddG4InTranscriptome(G4InTranscript,transcriptId, descriptionG4,informationsOfG4,localisationInTranscript, biotypeTranscript)
					G4InGenome=AddG4InGenome(G4InGenome, geneId, descriptionG4, localisationInTranscript)
					TranscriptPerG4=AddTranscriptPerG4(TranscriptPerG4, descriptionG4,transcriptId)

	ExtractionG4InTranscript(GITDIR+'/results/perChromosome', specie, chromosome, G4InTranscript)	
	ExtractionG4InGenome(GITDIR+'/results/perChromosome', specie, chromosome, G4InGenome)	
	ExtractionTranscriptPerG4(GITDIR+'/results/perChromosome', specie, chromosome, TranscriptPerG4)
	# ~ print "Done"
	
main()
