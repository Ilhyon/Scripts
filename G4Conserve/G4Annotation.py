#!/usr/bin/env python
# -*- coding: utf-8 -*-:

'''
Copyright Université de Sherbrooke, Département de Biochimie et Département d'Informatique

contact: sarah.belhamiti@usherbrooke.ca

This software is a computer program whose annote the G4 region of one chromosome of one specie.

---------------------------------------------------

``G4Annotation`` **module description**:

From a csv files which contain windows of G4 screener of a chromosome, this module annote
the G4 region for each transcript. The program create differents files output:

.. moduleauthor:: Sarah.Belhamiti

December 2017

'''

import sys
import os
import re
import recurentFunction as rF
from pprint import pprint
import argparse

def mean(liste):
   	return sum(liste)/len(liste)

def getChromosomalPositionForJunction(coord, strand, junLength,
									startFirstWindow,
									endFirstWindow):
	'''
		For more explanation abour those calculs go see 
		calculCoordsJunction.txt
	'''
	if strand == '1':
		coord = getChromosomalPositionForwardStrand(coord,
				junLength, startFirstWindow, endFirstWindow)
	elif strand == '-1':
		coord = getChromosomalPositionReverseStrand(coord,
				junLength, startFirstWindow, endFirstWindow)
	else :
		coord = None
	return coord

def isG4OnJunction(startG4, endG4, startBorder, endBorder):
	onJunction=False
	if ((startG4 < startBorder and endG4 > startBorder) or
		(startG4 > startBorder and endG4 < startBorder)):
		onJunction = True
	return onJunction

def getChromosomalPositionForwardStrand(position, EXTENSION,
									startIntron, endIntron):
	''' 
		Give real position (chromosomique) of start or end of a G4 
		around a junction. This function are for intron in the 
		forward strand.
	'''
	if position <= EXTENSION:
		# upstream squence (before the junction)
		position = startIntron - 1 - EXTENSION + position
	else:
		# downstream sequence (after the junction)
		position = endIntron + 1 - EXTENSION + position - 1
	return position

def getChromosomalPositionReverseStrand(position, EXTENSION,
									startIntron, endIntron):
	'''
		Give real position (chromosomique) of start or end of a G4 
		around a junction. This function are for intron in the 
		reverse strand.
	'''
	if position <= EXTENSION:
		# upstream squence (before the junction)
		position = startIntron + 1 + EXTENSION + 1 - position
	else:
		# downstream sequence (after the junction)
		position = endIntron + EXTENSION - position
	return position

def createIdG4(gene, startG4, endG4, strand):
	if strand == '1':
		return gene + '|' + str(startG4) + '|' + str(endG4) + '|' + strand
	elif strand == '-1':
		return gene + '|' + str(endG4) + '|' + str(startG4) + '|' + strand
	
def addWindowToG4Seq(g4Seq, windowSeq, step, windowLength):
	if len(windowSeq) < windowLength: 
		# if the sequence is shorter than the sequence of the windows
		# (generraly last windows above the thresolds)
		g4Seq += windowSeq[-(len(windowSeq)-(windowLength-step)):]
	else: #
		g4Seq += windowSeq[-step:] # take the stepsize added
	return g4Seq
	
def retrieveG4CoordByStrand(strand,
						  startWindow,
						  endWindow):
	if strand == '1':
		startG4 = startWindow
		endG4 = endWindow
	else :
		startG4 = endWindow
		endG4 = startWindow
	return startG4, endG4
	
def readLineG4Screener(line, StrandByGene, feature):
	'''
		There is many information in the output of G4RNA screener so 
		they will be contain in a dictionary for more readability.
		The last argument is the feature type (gene or junction).
		This allow to get 'good' names for keys in dico.
	'''
	words = line.split('\t') # get all element in a list
	dicoLine = {'Description' : words[1],
				'GeneID' : words[1].split('|')[0],
				'Strand' : StrandByGene.get(words[1].split('|')[0]),
				'cGcC' : float(words[2]),
				'g4H' : float(words[3]),
				'WindowSeq' : words[4],
				'WindowStart' : int(words[5]),
				'WindowEnd' : int(words[6]),
				'g4NN' : float(words[7])}
	if feature == 'Gene':
		dicoLine.update({feature+'Start' : int(words[1].split('|')[1]),
						feature+'End' : int(words[1].split('|')[2])})
	else:
		# for a junction, it is the coord of the intron
		dicoLine.update({feature+'Start' : int(words[1].split('|')[1].split('-')[0]),
						feature+'End' : int(words[1].split('|')[1].split('-')[1])})
	return dicoLine

def addCoordInlist(coordFeatureList, storageList):
	'''
		Browse the initial lsite of coords for a feature, parse them and
		then add them in a list.
		This list stores all coords of a transcript.
	'''
	for coord_feature in coordFeatureList:
		if coord_feature:
			startFeature = coord_feature.split('-')[0] 
			endFeature = coord_feature.split('-')[1]
			storageList.append(int(startFeature))
			storageList.append(int(endFeature))
	return storageList

def getCoordOfTranscript(start5 , end5 , start3 , end3,
							exonList, intronList, strand):
	''' 
		Get start and end of a transcript.
	'''
	listCoords = [] # store coords the transcript's features
	if not start5:
		# if the transcript doesn't have a 5' UTR
		if not end3:
			# if the transcript doesn't have a 3' UTR
			# coords are given by exons
			listCoords = addCoordInlist(exonList, listCoords)
			listCoords = addCoordInlist(intronList, listCoords)
			if strand == '1':
				# if forward gene then 5' < 3'
				inferiorCoord = min(listCoords)
				superiorCoord = max(listCoords)
			else:
				# if reverse gene then 5' > 3'
				inferiorCoord = max(listCoords)
				superiorCoord = min(listCoords)
		else: # if the transcript doesn't have 5' UTR
			# coords are given by exons and 3'
			listCoords = addCoordInlist(exonList, listCoords)
			listCoords = addCoordInlist(intronList, listCoords)
			if strand == '1':
				inferiorCoord = min(listCoords)
			else:
				inferiorCoord = max(listCoords)
			superiorCoord = int(end3)
	elif not start3: # if the transcript doesn't have a 3'UTR
		# border will be given by exons and 5'
		listCoords.append(start5)
		listCoords.append(end5)
		listCoords = addCoordInlist(exonList, listCoords)
		listCoords = addCoordInlist(intronList, listCoords)
		if strand == '1':
			superiorCoord = max(listCoords)
		else:
			superiorCoord = min(listCoords)
		inferiorCoord=int(start5)
	else: # if the transcript has a 5'UTR, a 3'UTR, exons and introns
		inferiorCoord = int(start5)
		superiorCoord = int(end3)
	return [inferiorCoord, superiorCoord]

def isG4InTranscript(strand, coordG4, coordsTranscript):	
	''' 
		Return a boolean : true if the pG4 is in a transcript or false
		on a the contriary.
	'''
	startG4 = int(coordG4[0])
	endG4 = int(coordG4[1])
	inferiorCoord = int(coordsTranscript[0])
	superiorCoord = int(coordsTranscript[1])
	inTranscript = False
	if strand == '1': # forward strand
		if (startG4 >= inferiorCoord and
			startG4 <= superiorCoord and
			endG4 >= inferiorCoord and
			endG4 <= superiorCoord):
			inTranscript = True
	else: # reverse strand
		if (startG4 <= inferiorCoord and
			startG4 >= superiorCoord and
			endG4 <= inferiorCoord and
			endG4 >= superiorCoord):
			inTranscript = True
	return inTranscript

def getLocalisationFeatureReverse(feature, featureList, localisation,
							startG4, endG4, biotype):
	for coord_feature in featureList:
			if coord_feature:
				# if transcript contain exon
				featureStart = int(coord_feature.split('-')[0])
				featureEnd = int(coord_feature.split('-')[1])
				if endG4 <= featureStart and endG4 >= featureEnd:
					# pG4 on an intron
					if startG4 >= featureEnd:
						# pG4 entierly in an intron
						localisation = biotype
						if feature == 'Intron' and biotype == 'CDS':
							localisation = 'Intron'
						elif feature == 'Intron and ' and biotype != 'CDS':
							localisation = 'IntronNC'
					else:
						# pG4 overlapping an exon and an intron
						if biotype == 'CDS':
							if feature == 'Exon':
								localisation = 'junction_'+biotype+'_Intron'
							else:
								localisation = 'junction_Intron_'+biotype
						else:
							if feature == 'Exon':
								localisation = 'junction_'+biotype+'_IntronNC'
							else:
								localisation = 'junction_IntronNC_'+biotype
	return localisation

def getLocalisationFeatureForward(feature, featureList, localisation,
							startG4, endG4, biotype):
	for coord_feature in featureList:
			if coord_feature:
				# if transcript contain exon
				featureStart = int(coord_feature.split('-')[0])
				featureEnd = int(coord_feature.split('-')[1])
				if startG4 >= featureStart and startG4 <= featureEnd:
					# pG4 on an intron
					if endG4 <= featureEnd:
						# pG4 entierly in an intron
						localisation = biotype
						if feature == 'Intron' and biotype == 'CDS':
							localisation = 'Intron'
						elif feature == 'Intron and ' and biotype != 'CDS':
							localisation = 'IntronNC'
					else:
						# pG4 overlapping an exon and an intron
						if biotype == 'CDS':
							if feature == 'Exon':
								localisation = 'junction_'+biotype+'_Intron'
							else:
								localisation = 'junction_Intron_'+biotype
						else:
							if feature == 'Exon':
								localisation = 'junction_'+biotype+'_IntronNC'
							else:
								localisation = 'junction_IntronNC_'+biotype
	return localisation

def getLocalisationForward(BiotypeByTranscript,
						ProteinCoding,
						transcriptId,
						coordG4, end5, start3,
						exonList, intronList):
	'''
		Return localisation of pG4 region in a transcript on forward
		strand.
	'''
	startG4 = int(coordG4[0]) 
	endG4 = int(coordG4[1])
	localisation = ''
	biotypeTranscript = BiotypeByTranscript.get(transcriptId)
	if biotypeTranscript in ProteinCoding:
		codant = 'CDS'
	else:
		codant = 'ExonNC'
	if end5 and startG4 <= int(end5): # if pG4 is positioned on the 5'UTR
		if endG4 <= int(end5): # if pG4 is enterly in the 5'UTR
			localisation='5'
		else: # if G4 is overlapping the 5'UTR
			localisation = 'junction_5_'+codant
	elif start3 and endG4 >= int(start3): # if pG4 is positioned on the 3'UTR
		if startG4 >= int(start3): # if pG4 is enterly in the 3'UTR
			localisation = '3'
		else: # if pG4 is overlapping the 3'UTR
			localisation = 'junction_'+codant+'_3' 
	else: # if pG4 is positioned in an intron
		localisation = getLocalisationFeatureForward('Exon', exonList,
				localisation, startG4, endG4, codant)
		localisation = getLocalisationFeatureForward('Intron', intronList,
				localisation, startG4, endG4, codant)
	return localisation

def getLocalisationReverse(BiotypeByTranscript,
						ProteinCoding,
						transcriptId,
						coordG4, end5, start3,
						exonList, intronList):
	'''
		Return localisation of pG4 region in a transcript on reverse
		strand.
	'''
	startG4 = int(coordG4[0])
	endG4 = int(coordG4[1])
	localisation = 'NAN'
	biotypeTranscript = BiotypeByTranscript.get(transcriptId)
	if biotypeTranscript in ProteinCoding:
		codant = 'CDS'
	else:
		codant = 'ExonNC'
	if end5 and endG4 >= int(end5): # if pG4 is positioned on the 5'UTR
		if startG4 >= int(end5): # if pG4 is enterly in the 5'UTR
			localisation = '5'
		else:	# if G4 is positioned in the 5' region with overlapping
			localisation = 'junction_5_'+codant
	elif start3 and startG4 <= int(start3): # if pG4 is positioned on the 3'UTR
		if endG4 <= int(start3): # if pG4 is enterly in the 3'UTR
			localisation = '3'
		else: # if pG4 is overlapping the 3'UTR
			localisation = 'junction_'+codant+'_3' 
	else: # if pG4 is positioned in an intron
		localisation = getLocalisationFeatureReverse('Exon', exonList,
				localisation, startG4, endG4, codant)
		localisation = getLocalisationFeatureReverse('Intron', intronList,
				localisation, startG4, endG4, codant)
	return localisation



def getLocalisationG4InJunction(BiotypeByTranscript,
								ProteinCoding,
								transcriptId):
	localisation='NAJ'
	biotypeTranscript = BiotypeByTranscript.get(transcriptId)
	if biotypeTranscript in ProteinCoding:
		codant='CDS'
	else:
		codant='ExonNC'
	return 'junction_'+codant+'_'+codant
	
def getLocalisationG4InTranscript(BiotypeByTranscript,
								ProteinCoding,
								transcriptId,
								coordG4, end5, start3,
								exonList, intronList, strand):
	''' 
		Depending on the strand, the localisation of the pG4
		is found.
	'''
	if strand == '1': # if gene on a forward strand
		localisation = getLocalisationForward(BiotypeByTranscript,
											ProteinCoding,
											transcriptId,
											coordG4, end5, start3,
											exonList, intronList)
	else: # if gene on a reverse strand
		localisation = getLocalisationReverse(BiotypeByTranscript,
											ProteinCoding,
											transcriptId,
											coordG4, end5, start3,
											exonList, intronList)
	return localisation

def addG4InTranscriptome(G4InTranscript, transcriptId,
						descriptionG4, informationsOfG4,
						localisationInTranscript, biotypeTranscript):
	''' 
		This fonction add for each pG4 the informations from 
		G4RNA screener (cGcC, g4H, sequence and G4NN) to the 
		localisation of the pG4 and the biotype of the transcript.
		Those Information are stored in a dico :
		{idG4 : informations}
	'''
	value = []
	for information in informationsOfG4:
		value.append(information)
	value.append(localisationInTranscript)
	value.append(biotypeTranscript)
	headerG4 = transcriptId+'|'+descriptionG4 # uniq header pG4
	if headerG4 not in G4InTranscript:
		G4InTranscript[headerG4] = value
	return G4InTranscript

def addG4InGenome(G4InGenome, geneId, descriptionG4, localisationInTranscript):
	''' 
		This fonction add for each gene all pG4 in it and 
		all localisation possible.
		{idG4 : info}
	'''
	headerG4 = geneId+':'+descriptionG4
	if headerG4 not in G4InGenome:
		G4InGenome[headerG4] = localisationInTranscript
	else:
		localisationInGene=G4InGenome.get(headerG4)
		if localisationInGene :
			if localisationInTranscript not in localisationInGene:
				G4InGenome[headerG4] = localisationInGene+';'+localisationInTranscript
	return G4InGenome

def addTranscriptPerG4(TranscriptPerG4, descriptionG4,transcriptId): 
	''' 
		This fonction add for each pG4 all transcript that contain it.
		{idG4 : info}
	'''
	if descriptionG4 not in TranscriptPerG4: 
		TranscriptPerG4[descriptionG4] = transcriptId
	else:
		if TranscriptPerG4:
			lastTranscript = TranscriptPerG4.get(descriptionG4)
			TranscriptPerG4[descriptionG4] = lastTranscript +'-'+transcriptId
	return TranscriptPerG4

def isJuncionInTranscript(coordG4,intronList):
	''' 
		Return a boolean : true if a junction with a pG4 is located 
		in a transcript and false on the contriary.
	'''
	startG4 = coordG4[0]
	endG4 = coordG4[1]
	inTranscript = False
	for couple_Intron in intronList:
		if couple_Intron:
			endIntron = int(couple_Intron.split('-')[1])
			if startG4 <= endIntron and endIntron<= endG4:
				inTranscript = True
	return inTranscript

def getUTRCoordByStrand(words, dicoLine):
	if dicoLine['Strand'] == '1':
		dicoLine.update({'start5' : words[7],
						'end5' : words[8],
						'start3' : words[9],
						'end3' : words[10].rstrip()})
	else:
		dicoLine.update({'start5' : words[8],
						'end5' : words[7],
						'start3' : words[10].rstrip(),
						'end3' : words[9]})
	return dicoLine

def readLineIndex(line):
	words = line.split('|')
	dicoLine = {'idTr' : words[0],
				'idGene' : words[1],
				'Chromosome' : words[2],
				'Strand' : words[3],
				'geneBiotype' : words[4],
				'exonList' : words[5].split(';'),
				'intronList' : words[6].split(';')}
	dicoLine.update(getUTRCoordByStrand(words, dicoLine))
	return dicoLine

def extractionTranscriptPerG4(directory, specie, chromosome, TranscriptPerG4):
	output= open(directory+'/'+specie+'_chr'+chromosome+'_TranscriptPerG4.txt','w')
	output.write('CoordonnéesG4\tTranscript(s)\n')
	for key,value in TranscriptPerG4.items():
		output.write(key+'\t'+value+'\n')
	output.close()
	
def extractionG4InGenome(directory, specie, chromosome, G4InGenome):
	output= open(directory+'/'+specie+'_chr'+chromosome+'_G4InGenome.txt','w')
	output.write('InfoG4ByGene\tLocalisation(s)\n')
	for key,value in G4InGenome.items():
		output.write(key+'\t'+value+'\n')
	output.close()
	
def extractionG4InTranscript(directory, specie, chromosome, G4InTranscript):
	output= open(directory+'/'+specie+'_chr'+chromosome+'_G4InTranscript.txt','w')
	output.write('InfoG4ByTranscript\tcGcC\tG4Hunter\tsequenceG4\tG4NN\tlocalisation\ttranscriptBiotype\n')
	for key,value in G4InTranscript.items():
		if None not in value:
			# because some transcriptID from ensembl don't 
			# contains info of biotype (as ENST00000604369)
			output.write(key+'\t'+'\t'.join(value)+'\n')
	output.close()

def getInfoAboutpG4(index, BiotypeByTranscript, listeG4InGeneEntire,
					G4InTranscript, G4InGenome, TranscriptPerG4,
					AnnotationTranscript, G4DetectedInGene, ProteinCoding,
					feature):
	inputfile = open(index,'r') # file opening for reading
	for line in inputfile: # for each transcript
		l = readLineIndex(line)
		biotypeTranscript = BiotypeByTranscript.get(l['idTr'])
		borderTranscript = getCoordOfTranscript(l['start5'], l['end5'],
												l['start3'], l['end3'],
												l['exonList'],
												l['intronList'],
												l['Strand'])
		if l['idGene'] in listeG4InGeneEntire:
			listeG4InGene = listeG4InGeneEntire[ l['idGene'] ]
			for G4InGene in listeG4InGene:
				startG4 = int(G4InGene.split('|')[1])
				endG4 = int(G4InGene.split('|')[2])
				coordG4 = [startG4,endG4]
				if feature == 'Transcript':
					g4inTranscript = isG4InTranscript(l['Strand'],
									coordG4, borderTranscript)
				elif feature == 'Junction':
					g4inTranscript = isJuncionInTranscript(coordG4,
											l['intronList'])
				annotationTranscript = AnnotationTranscript[ l['idTr'] ]
				if g4inTranscript and annotationTranscript: # PBM ICI COORD
					informationsOfG4 = list(G4DetectedInGene[G4InGene])
					if feature == 'Transcript':
						localisationInTranscript = getLocalisationG4InTranscript(BiotypeByTranscript,
											ProteinCoding, l['idTr'],
											coordG4, l['end5'], l['start3'],
											l['exonList'], l['intronList'],
											l['Strand'])
					elif feature == 'Junction':
						localisationInTranscript = getLocalisationG4InJunction(BiotypeByTranscript,
											ProteinCoding, l['idTr'])
					descriptionG4 = l['Chromosome']+':'+str(startG4)+'-'+str(endG4)+'|'+l['Strand']
					G4InTranscript = addG4InTranscriptome(G4InTranscript,
									l['idTr'],
									descriptionG4,
									informationsOfG4,
									localisationInTranscript,
									biotypeTranscript)
					G4InGenome = addG4InGenome(G4InGenome,
												l['idGene'],
												descriptionG4,
												localisationInTranscript)
					TranscriptPerG4 = addTranscriptPerG4(TranscriptPerG4, descriptionG4,l['idTr'])
	inputfile.close()
	return G4InTranscript, G4InGenome, TranscriptPerG4

def getlisteG4InGene(G4Detected):
	'''
		Return a dictionary contain all pG4 of a gene.
		{idGene : pG4 list}
	'''
	listeG4InGene = {}
	for descriptionG4 in G4Detected:
		geneID = descriptionG4.split('|')[0]
		if geneID not in listeG4InGene:	
			listeG4 = []
		else:
			listeG4 = listeG4InGene.get(geneID)
		listeG4.append(descriptionG4)
		listeG4InGene[geneID] = listeG4
	return listeG4InGene

def returnG4InGene(G4DetectedInGene,
					inputfile,
					dicoParam,
					StrandByGene):
	'''
		Returns a dictionnary with informations of predicted region G4.
		Those informations contains : score cgCc, scote G4Hunter,
		sequence and score G4NN.
		This function browse all windows returned by G4RNA Screener and 
		keep only those over the thresholds. It also merge overlapping
		windows.
	'''
	oldPassed = False # boolean about the previous windows
	passed = False # boolean : true if over threshold, elsewise false
	inputfile = open(inputfile,'r')
	for line in inputfile:
		if (re.search('^[0-9]', line)):
			# if the line is not the header of the file 
			dicoLine = readLineG4Screener(line, StrandByGene, 'Gene')
			if (dicoLine['cGcC'] >= dicoParam['cGcC'] and
				dicoLine['g4H'] >= dicoParam['g4H'] and
				dicoLine['g4NN'] >= dicoParam['g4NN']):
				# if window over thresholds
				passed = True # we are over thresholds
				if (oldPassed != passed):
					# first windows over the threshold
					descriptionOverThreshold = dicoLine['Description'] 
					# assignation of description for this window
					sequenceG4 = dicoLine['WindowSeq']
					oldPassed = passed # update
					if dicoLine['Strand'] == '1': # if forward gene
						startG4 = dicoLine['WindowStart']
						endG4 = dicoLine['WindowEnd']
					else: 
						# cf calcul 1 from calculCoords
						startG4 = dicoLine['GeneEnd'] - \
								(dicoLine['WindowStart'] -\
								dicoLine['GeneStart'])
						endG4 = startG4 - (len(dicoLine['WindowSeq'])) + 1
					# initialization of a list for eahc score
					listeCGcC = [ dicoLine['cGcC'] ]
					listeG4Hunter = [ dicoLine['g4H'] ]
					listeG4NN = [ dicoLine['g4NN'] ]
				else:
					# not the first windows above the thresholds
					sequenceG4 = addWindowToG4Seq(sequenceG4,
								dicoLine['WindowSeq'],
								dicoParam['Step'], dicoParam['Window'])
					if dicoLine['Strand'] == '1' :
						endG4 = dicoLine['WindowEnd']
					else :
						endG4 = startG4 - (len(sequenceG4)) + 1
					listeCGcC.append(dicoLine['cGcC'])
					listeG4Hunter.append(dicoLine['g4H'])
					listeG4NN.append(dicoLine['g4NN'])
			if (dicoLine['cGcC'] < dicoParam['cGcC'] or
				dicoLine['g4H'] < dicoParam['g4H'] or
				dicoLine['g4NN'] < dicoParam['g4NN'] or
				descriptionOverThreshold != dicoLine['Description']):
				# one of the score is under his threshold 
				# or this windows is from another gene
				passed=False # update
				if (oldPassed != passed ):
					# last windows before under the thresolds
					meanCGcC = mean(listeCGcC) 
					meanG4Hunter = mean(listeG4Hunter)
					meanG4NN = mean(listeG4NN)
					oldPassed = passed # update
					headerG4 = createIdG4(dicoLine['GeneID'],
								startG4, endG4, dicoLine['Strand'])
					if headerG4 not in G4DetectedInGene and dicoLine['Strand']: 
						G4DetectedInGene[headerG4] = [str(meanCGcC),
													str(meanG4Hunter),
													sequenceG4,
													str(meanG4NN)]	
	inputfile.close()
	return G4DetectedInGene

def returnG4InJunction(G4DetectedInJunction,
						inputfile,
						dicoParam,
						StrandByGene):
	'''
		Returns a dictionnary with informations of predicted region G4.
		Those informations contains : score cgCc, scote G4Hunter,
		sequence and score G4NN.
		This function browse all windows returned by G4RNA Screener and 
		keep only those over the thresholds. It also merge overlapping
		windows.
	'''
	oldPassed = False # boolean about the previous windows
	passed = False # boolean : true if over threshold, elsewise false
	descriptionOverThreshold = ''
	inputfile = open(inputfile,'r')	# file opening for reading
	for line in inputfile:	# for each file csv create by G4 screener
		if (re.search('^[0-9]', line)): # if the line is not the header of the file 
			dicoLine = readLineG4Screener(line, StrandByGene, 'Junction')
			if (dicoLine['cGcC'] >= dicoParam['cGcC'] and
				dicoLine['g4H'] >= dicoParam['g4H'] and
				dicoLine['g4NN'] >= dicoParam['g4NN'] and 
				dicoLine['WindowStart']):
				# window over thresholds
				onJunction = False
				passed = True
				if (oldPassed != passed):
					# first windows over thresholds
					descriptionOverThreshold = dicoLine['Description'] 
					# assignation of description for this window
					sequenceG4 = dicoLine['WindowSeq']
					oldPassed = passed # update
					startFirstWindow = dicoLine['JunctionStart']
					endFirstWindow = dicoLine['JunctionEnd']
					startG4 = dicoLine['WindowStart']
					endG4 = dicoLine['WindowEnd']
					# the pG4 start and end like the window
					# initialization of a list for eahc score
					listeCGcC = [ dicoLine['cGcC'] ]
					listeG4Hunter = [ dicoLine['g4H'] ]
					listeG4NN = [ dicoLine['g4NN'] ]
				else:
					# not the first windows above the thresholds
					sequenceG4 = addWindowToG4Seq(sequenceG4,
								dicoLine['WindowSeq'],
								dicoParam['Step'],
								dicoParam['Window'])
					endG4 = dicoLine['WindowEnd']
					listeCGcC.append(dicoLine['cGcC'])
					listeG4Hunter.append(dicoLine['g4H'])
					listeG4NN.append(dicoLine['g4NN'])
			if (dicoLine['cGcC'] < dicoParam['cGcC'] or
				dicoLine['g4H'] < dicoParam['g4H'] or
				dicoLine['g4NN'] < dicoParam['g4NN'] or
				descriptionOverThreshold != dicoLine['Description']):
				# one of the score is under his threshold 
				# or this windows is from another gene
				passed = False
				if (oldPassed != passed ):
					# last windows before under the thresolds
					meanCGcC = mean(listeCGcC) 
					meanG4Hunter = mean(listeG4Hunter)
					meanG4NN = mean(listeG4NN)
					oldPassed=passed # update
					startG4 = getChromosomalPositionForJunction(startG4,
									dicoLine['Strand'], dicoParam['Window'],
									startFirstWindow, endFirstWindow)
					endG4 = getChromosomalPositionForJunction(endG4,
									dicoLine['Strand'], dicoParam['Window'],
									startFirstWindow, endFirstWindow)
					headerG4 = createIdG4(dicoLine['GeneID'],
								startG4, endG4, dicoLine['Strand'])
					onJonction = isG4OnJunction(startG4, endG4,
								dicoLine['JunctionStart'], dicoLine['JunctionEnd'])	
					if headerG4 not in G4DetectedInJunction and onJonction:
						G4DetectedInJunction[headerG4] = [str(meanCGcC),
														str(meanG4Hunter),
														sequenceG4,
														str(meanG4NN)]
	inputfile.close()
	return G4DetectedInJunction

def createDictionaryStrandByGene(filename): 
	''' 
		Create dictionary with the strand of all gene.
		Dico : {gene's id : gene's strand}
	'''
	dico = {}
	inputfile = open(filename,'r')
	for line in inputfile:
		words = line.split('|')
		gene = words[1]
		strand = words[3]
		if gene not in dico:
			dico[gene]=strand
	inputfile.close()
	return dico

def createListCodingProtein():
	codingProtein=['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene',
					'IG_M_gene', 'IG_V_gene', 'IG_Z_gene',
					'nonsense_mediated_decay', 'protein_coding',
					'nontranslating_CDS', 'non_stop_decay',
					'TR_C_gene', 'TR_D_gene', 'TR_gene', 
					'TR_J_gene', 'TR_V_gene']
	return codingProtein

def createDicoParam(arg):
	dicoParam = {'g4H' : float(arg.THRESHOLD_G4H),
				'cGcC' : float(arg.THRESHOLD_CGCC),
				'g4NN' : float(arg.THRESHOLD_G4NN),
				'Extension' : int(arg.EXTENSION),
				'Window' : int(arg.EXTENSION),
				'Step' : int(arg.STEP)}
	return dicoParam

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	GITDIR=os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR+'data')
	parser.add_argument ('-chr', '--chromosome', default = 'X')
	parser.add_argument ('-specie', '--specie', default = 'HS')
	parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-E', '--EXTENSION', default = 100)
	parser.add_argument ('-W', '--WINDOW', default = 60)
	parser.add_argument ('-S', '--STEP', default = 10)
	return parser

def main () :
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path # directory containing all chromosome
	chromosome = arg.chromosome # chromosome to analyze
	specie = arg.specie # specie to analyse
	dicoParam = createDicoParam(arg)
	GITDIR = os.getcwd()+'/'
	ProteinCoding = createListCodingProtein()
	G4DetectedInGene = {}
	G4DetectedInJunction = {}
	G4InTranscript = {}
	G4InGenome = {}
	TranscriptPerG4 = {}
	directory = path+'/chr'+chromosome
	# variable directory which contain the data for this chromosome
	index = directory+'/'+specie+'_transcript_unspliced_chr'+chromosome+'_Index.txt'
	# file which contain info by transcript for this chromosome
	indexBiotypeTranscript = path+'/transcriptType/transcriptType_chr'+chromosome
	# file which contain biotype of transcript for this chromosome
	print "Chromosome " + chromosome
	BiotypeByTranscript = rF.createDictionaryBiotypeByTranscript(indexBiotypeTranscript)
	StrandByGene = createDictionaryStrandByGene(index)
	AnnotationTranscript = rF.GetAnnotationTranscript(index, 
							ProteinCoding, BiotypeByTranscript)
	# get g4 from the ouput of G4RNA Screener
	for path, dirs, files in os.walk(directory):
		# for each element of the directory to passed
		for filename in files: # for each files (all csv)
			inputfile = directory+'/'+filename
			if ('gene_unspliced' in filename and '.csv' in filename ):	
				G4DetectedInGene = returnG4InGene(G4DetectedInGene, 
									inputfile, dicoParam,
									StrandByGene)
			elif ('transcript_unspliced' in filename and '.csv' in filename):
				G4DetectedInJunction = returnG4InJunction(G4DetectedInJunction,
										inputfile, dicoParam,
										StrandByGene)
	listeG4InGeneEntire = getlisteG4InGene(G4DetectedInGene)
	listeG4InGeneJunction = getlisteG4InGene(G4DetectedInJunction)
	G4InTranscript, G4InGenome, TranscriptPerG4 = getInfoAboutpG4(index,
					BiotypeByTranscript, listeG4InGeneEntire,
					G4InTranscript, G4InGenome, TranscriptPerG4,
					AnnotationTranscript, G4DetectedInGene, ProteinCoding,
					'Transcript')
	G4InTranscript, G4InGenome, TranscriptPerG4 = getInfoAboutpG4(index,
					BiotypeByTranscript, listeG4InGeneJunction,
					G4InTranscript, G4InGenome, TranscriptPerG4,
					AnnotationTranscript, G4DetectedInJunction, ProteinCoding,
					'Junction')
	extractionG4InTranscript(GITDIR+'results/perChromosome', specie, chromosome, G4InTranscript)	
	extractionG4InGenome(GITDIR+'results/perChromosome', specie, chromosome, G4InGenome)	
	extractionTranscriptPerG4(GITDIR+'results/perChromosome', specie, chromosome, TranscriptPerG4)
	print "\t Done."
	
main()
