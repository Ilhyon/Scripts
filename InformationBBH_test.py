#----------------------------------------------------------------------#
def GetCommonInformation(localisationsSp1, localisationsSp2, biotypesSp1, biotypesSp2):
	""" Get the commons informations between two list for the biotypes and the localisations,
		if there is no commons informations, the list will stay empty
		Parameters
	    ----------
	    localisationsSp1 :	list of localisations from the specie 1
		localisationsSp2 :	list of localisations from the specie 2
		biotypesSp1 : list of biotypes from the specie 1
		biotypesSp2 : list of biotypes from the specie 2
		
	    Returns
	    -------
	    commonLocalisations : list of commons localisations of a G4 in a gene/transcript between the 2 specie
	    commonBiotypes : list of commons biotypes of a G4 in a gene/transcript between the 2 specie
	"""
	commonLocalisations = []
	commonBiotypes = []
	if list(set(localisationsSp1).intersection(localisationsSp2)) :
		commonLocalisations = list(set(localisationsSp1).intersection(localisationsSp2))
	if list(set(biotypesSp1).intersection(biotypesSp2)) :
		commonBiotypes = list(set(biotypesSp1).intersection(biotypesSp2))
	return(commonLocalisations, commonBiotypes)
#----------------------------------------------------------------------#
def GetLocalisationsAndBiotype(element1, element2, Dico):
	""" Get the localisations and the biotypes in a dictionary of dictionary 
		thanks to 2 elements
		Parameters
	    ----------
	    element1 :	string, element that bust be present in the first level of the dictionary
		element2 :	string, element that bust be present in the second level of the dictionary
		Dico :	{idTranscript/Gene : {idG4 : {"Localisation" : [localisation]
											  "Biotype" : [Biotype]}}}
	    Returns
	    -------
	    localisations : list of localisations of a G4 in a gene/transcript
	    biotypes : list of biotypes of a G4 in a gene/transcript
	"""
	if element1 in Dico :
		if element2 in Dico[element1] :
			localisations = Dico[element1][element2]["Localisations"]
			biotypes = Dico[element1][element2]["Biotype"]
		else :
			print "Patate 4"
	else :
		print "Patate 5"
	return(localisations, biotypes)
#----------------------------------------------------------------------#
def AddInformationsOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleOrtho, commonLocalisations, commonBiotypes, typeInfo):
	""" Adds the informations to the final dictionary containing All 
		infromation about the BBH if there is some orthology 
		Parameters
	    ----------
	    DicoInfoBBHLevel :
		idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		coupleOrtho : string, id of the couple like geneSp1|geneSp2 or transcriptSp1|transcriptSp2
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {"Orthology" : {coupleOrtho : {"Common Localisations" : [commonLocalisations]
																	   "Common Biotypes" : commonBiotypes}}}}
	"""
	if idBBH not in DicoInfoBBHLevel :
		# First time the couple of BBH is encounter
		# so we create all information relative to it 
		DicoInfoBBHLevel[idBBH] = {typeInfo : {}}
		DicoInfoBBHLevel[idBBH][coupleOrtho] = {}
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and coupleOrtho not in DicoInfoBBHLevel[idBBH][typeInfo] :
		# the couple of BBH is already in the dictionary but 
		# not the couple there homologues's levelS
		DicoInfoBBHLevel[idBBH][coupleOrtho] = {}
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and coupleOrtho in DicoInfoBBHLevel[idBBH] :
		# the couple of BBH is already in the dictionary and also the couple of homologue
		# yet we still have a chance to get some new informations so we need to check it
		if list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)) :
			# if there is new localisations we add them
			DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Localisations"].append(list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)))
		if list(set(DicoInfoBBHLevel[idBBH]["Common Biotypes"]) ^ set(commonBiotypes)) :
			# if there is new biotypes we add them
			DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Biotypes"].append(commonBiotypes)
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def AddInfoNoOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, newList, commonLocalisations, commonBiotypes, typeInfo):
	""" Adds the informations to the final dictionary containing All infromation 
		about the BBH if there is no orthology links
		Parameters
	    ----------
	    DicoInfoBBHLevel :
		idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		idNoOrthology : list of string, id of the gene that are not orthologues
						but they still have common informations in common
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {"Orthology" : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
																	   "Common Biotypes" : commonBiotypes}}}}
	"""
	if idBBH not in DicoInfoBBHLevel :
		# First time the couple of BBH is encounter
		# so we create all information relative to it 
		DicoInfoBBHLevel[idBBH] = {typeInfo : {}}
		DicoInfoBBHLevel[idBBH][typeInfo][newList] = {}
		DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and newList not in DicoInfoBBHLevel[idBBH][typeInfo] :
		cpt = 0
		# the couple of BBH is already in the dictionary but 
		# not the id that are not orthologues, this can be du to many reasons
		# maybe there just one id in plus the list or maybe the list doesn't exist
		for existingList in DicoInfoBBHLevel[idBBH][typeInfo]:
			if set(existingList).issubset(set(newList)) :
				cpt +=1
				# the existingList is a subset of newList
				newComLoca = list(set(commonLocalisations) ^ (DicoInfoBBHLevel[idBBH][typeInfo][existingList]["Common Localisations"]))
				newComBiot = list(set(commonBiotypes) ^ (DicoInfoBBHLevel[idBBH][typeInfo][existingList]["Common Biotypes"]))
				DicoInfoBBHLevel[idBBH][typeInfo][new_key] = dictionary.pop(DicoInfoBBHLevel[idBBH][typeInfo][existingList])
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"].append(newComLoca)
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"].append(newComBiot)
			else :
				# the list of id is new
				DicoInfoBBHLevel[idBBH][typeInfo][newList] = {}
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"] = commonLocalisations
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"] = commonBiotypes
		if cpt > 1 :
			print "Patate 7"		
	elif idBBH in DicoInfoBBHLevel and newList in DicoInfoBBHLevel[idBBH][typeInfo] :
		# the couple of BBH is already in the dictionary and also the couple of homologue
		# yet we still have a chance to get some new informations so we need to check it
		if list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)) :
			# if there is new localisations we add them
			DicoInfoBBHLevel[idBBH][idNoOrthology]["Common Localisations"].append(list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)))
		if list(set(DicoInfoBBHLevel[idBBH]["Common Biotypes"]) ^ set(commonBiotypes)) :
			# if there is new biotypes we add them
			DicoInfoBBHLevel[idBBH][idNoOrthology]["Common Biotypes"].append(commonBiotypes)
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def AddingInformations(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleHomologue, commonLocalisations, commonBiotypes, typeInfo):
	""" Adds the informations to the final dictionary containing All infromation about the BBH
		Parameters
	    ----------
	    DicoInfoBBHLevel :
		idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		coupleHomologue : 
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
	"""
	idBBH = str(idG4Sp1+"|"+idG4Sp2) # creation of the id of the couple of BBH
	
	if typeInfo == "Orthology" :
		DicoInfoBBHLevel = AddInformationsOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleHomologue, commonLocalisations, commonBiotypes, typeInfo)
	else if typeInfo == "No orthology" :
		
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetInfoNoOrthology(levelSp1, nonHomologueslevelSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel):
	""" Here we test if there is common information between gene that hasn't orthology links
		Parameters
	    ----------
	    nonHomologueslevelSp2 : list of string, list of the gene that are not orthologs to the gene of 
								the specie 1 where there is the G4 BBH we are looking at for
		levelSp1 : string, transcript/gene from the specie 1
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		level : string, to choose if we want the gene level or the transcript level
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleNoOrth : {"Common Localisations" : [commonLocalisations]
													 "Common Biotypes" : commonBiotypes}}}
	"""
	idNoOrthology = [levelSp1]
	commonLocalisationsNoOrth = []
	commonBiotypesNoOrth[]
	for levelSp2 in nonHomologueslevelSp2 : 
		# retrivial of the localisations and biotype of the two specie
		# from the gene/transcript that are homologues
		localisationsSp1WithoutH, biotypesSp1WithoutH = GetLocalisationsAndBiotype(idG4Sp1, levelSp1, DicoInfoG4Sp1Level)
		localisationsSp2WithoutH, biotypesSp2WithoutH = GetLocalisationsAndBiotype(idG4Sp2, levelSp2, DicoInfoG4Sp2Level)
		commonLocalisationsWithoutH, commonBiotypesWithoutH = GetCommonInformation(localisationsSp1WithoutH, localisationsSp1WithoutH, biotypesSp1WithoutH, biotypesSp1WithoutH)
		if commonLocalisationsWithoutH or commonBiotypesWithoutH :
			commonLevelSp2.append(levelSp2)
			commonLocalisationsNoOrth.append(commonLocalisationsWithoutH)
			commonBiotypesNoOrth.append(commonBiotypesWithoutH)
	commonLocalisationsNoOrth = list(set(commonLocalisationsNoOrth))
	commonBiotypesNoOrth =  list(set(commonBiotypesNoOrth))
	DicoInfoBBHLevel = AddingInformations(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, idNoOrthology, commonLocalisationsNoOrth, commonBiotypesNoOrth, "No orthology")
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetInfoHomology(HlevelSp1, levelSp1, levelsSp2, idG4Sp2, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel):
	""" Here we test many conditions (with or without homology), then we 
		retrieve the informations relative to them
		Parameters
	    ----------
	    HlevelSp1 : list of string, list of the gene that are homologue to the gene of 
					the specie 1 where there is the G4 BBH we are looking at for
		levelSp1 : string, transcript/gene from the specie 1
		levelsSp2 : list of string, list of transcript/gene from the specie 2
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		level : string, to choose if we want the gene level or the transcript level
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
	"""
	if list(set(HlevelSp1).intersection(levelsSp2)) :
		# if there an intersection between the homologues of the gene/transcript of specie 1
		# and with the levelS of the specie 2 where the G4 is
		
		for levelSp2 in list(set(HlevelSp1).intersection(levelsSp2)) :
			# retrivial of the localisations and biotype of the two specie
			# from the gene/transcript that are homologues
			localisationsSp1WithH, biotypesSp1WithH = GetLocalisationsAndBiotype(idG4Sp1, levelSp1, DicoInfoG4Sp1Level)
			localisationsSp2WithH, biotypesSp2WithH = GetLocalisationsAndBiotype(idG4Sp2, levelSp2, DicoInfoG4Sp2Level)
			commonLocalisationsWithH, commonBiotypesWithH = GetCommonInformation(localisationsSp1WithH, localisationsSp2WithH, biotypesSp1WithH, biotypesSp2WithH)
			coupleOrthologue = str(levelSp1+"|"+levelSp2)
			DicoInfoBBHLevel = AddingInformations(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleOrthologue, commonLocalisationsWithH, commonBiotypesWithH, "Orthology")
			
		nonHomologueslevelSp2 = list(set(HlevelSp1) ^ set(levelsSp2))
		if nonHomologueslevelSp2 : 
		# retrivial of the list of G4 that are not homologues (unique in the two lists)
			DicoInfoBBHLevel = GetInfoNoOrthology(levelSp1, nonHomologueslevelSp2, idG4Sp2, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel)
			
	else :
		# there is no element in comon => no orthologues between the gene/transcript
		# where the G4 BBH are
		DicoInfoBBHLevel = GetInfoNoOrthology(levelSp1, levelsSp2, idG4Sp2, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel)
return(DicoInfoBBHLevel)
