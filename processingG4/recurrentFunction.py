#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

def setUpperLetter(words):
    """Get all first letter of a name and set them in upper case.

    .. admonition:: Example of changement.

		The specie escherichia_coli_str_k_12_substr_mg1655 gave ECSK1SM

    :param words: words that need to get their first  letter in upper case,
	:type words: string

    :returns: ini, initials letter upper cased.
    :rtype: string
    """
    words = words.split("_")
    letters = [word[0] for word in words]
    ini = "".join(letters)
    ini = ini.upper()
    return ini

def getShortSpName(sp):
    sp = sp.split("_")
    firstLetter = sp[0][0]
    firstLetter = firstLetter.upper()
    other = sp[1][0:3]
    shortName = "".join([firstLetter, other])
    return shortName

def createDicoParam(arg):
	"""Retrieves arguments and put them in a dictionary.

	:param arg: contains all arguments given to the script, those are principaly
		parameters from G4RNA Screener.
	:type arg: arg_parser

	:returns: dicoParam, contains all arguments given to the script.
	:rtype: dictionary
	"""
	dicoParam = {"G4H" : float(arg.THRESHOLD_G4H),
				"cGcC" : float(arg.THRESHOLD_CGCC),
				"G4NN" : float(arg.THRESHOLD_G4NN),
				"junctionLength" : int(arg.EXTENSION),
				"windowLength" : int(arg.EXTENSION),
				"step" : int(arg.STEP)}
	return dicoParam

def createScoreDico():
    dico = {'G4NN': {'G4NN' : 0.0,
                    'cGcC' : 4.5,
                    'G4H' : 0.9,},
            'G4H' : {'G4NN' : 0.5,
                    'cGcC' : 4.5,
                    'G4H' : -4.6},
            'cGcC' : {'G4NN' : 0.5,
                    'cGcC' : 0.0,
                    'G4H' : 0.9,}}
    return dico

def createDicoVenn():
    dico = {'G4NN' : {'g4NN' : 0.5,
                        'cGcC' : -4.5,
                        'g4H' : -4.5},
            'G4NN-G4H' : {'g4NN' : 0.5,
                            'cGcC' : -4.5,
                            'g4H' : 0.9},
            'G4NN-cGcC' : {'g4NN' : 0.5,
                            'cGcC' : 4.5,
                            'g4H' : -4.5},
            'G4H' : {'g4NN' : 0,
                    'cGcC' : -4.5,
                    'g4H' : 0.9},
            'G4H-cGcC' : {'g4NN' : 0,
                        'cGcC' : 4.5,
                        'g4H' : 0.9},
            'cGcC' : {'g4NN' : 0,
                    'cGcC' : 4.5,
                    'g4H' : -4.5},
            'cGcC-G4H-G4nn' : {'g4NN' : 0.5,
                                'cGcC' : 4.5,
                                'g4H' : 0.9}}
    return dico
