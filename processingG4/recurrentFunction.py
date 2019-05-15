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

def createDicoParam(arg):
	"""Retrieves arguments and put them in a dictionary.

	:param arg: contains all arguments given to the script, those are principaly
		parameters from G4RNA Screener.
	:type arg: arg_parser

	:returns: dicoParam, contains all arguments given to the script.
	:rtype: dictionary
	"""
	dicoParam = {"g4H" : float(arg.THRESHOLD_G4H),
				"cGcC" : float(arg.THRESHOLD_CGCC),
				"g4NN" : float(arg.THRESHOLD_G4NN),
				"junctionLength" : int(arg.EXTENSION),
				"windowLength" : int(arg.EXTENSION),
				"step" : int(arg.STEP)}
	return dicoParam
