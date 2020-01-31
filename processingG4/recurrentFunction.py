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
				"junctionLength" : int(arg.junctionLength),
				"windowLength" : int(arg.junctionLength),
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

def createDicoSp():
    dico = {'Aaeo' : 1,
            'Mtub' : 2,
            'Tthe' : 3,
            'Caur' : 4,
            'Mpne' : 5,
            'Saur' : 6,
            'Bsub' : 7,
            'Efae' : 8,
            'Spne' : 9,
            'Ctra' : 10,
            'Bbur' : 11,
            'Cjej' : 12,
            'Mxan' : 13,
            'Gsul' : 14,
            'Babo' : 15,
            'Wend' : 16,
            'Apha' : 17,
            'Nmen' : 18,
            'Lpne' : 19,
            'Ftul' : 20,
            'Vcho' : 21,
            'Hinf' : 22,
            'Ypes' : 23,
            'Ecol' : 24,
            'Csym' : 25,
            'Ckor' : 26,
            'Nequ' : 27,
            'Paer' : 28,
            'Ssol' : 29,
            'Taci' : 30,
            'Phor' : 31,
            'Msmi' : 32,
            'Aful' : 33,
            'Hbut' : 34,
            'Hsal' : 35,
            'Mace' : 36,
            'Lmaj' : 37,
            'Crei' : 38,
            'Ppat' : 39,
            'Osat' : 40,
            'Slyc' : 41,
            'Vvin' : 42,
            'Atha' : 43,
            'Ddis' : 44,
            'Spom' : 45,
            'Scer' : 46,
            'Anid' : 47,
            'Ncra' : 48,
            'Cele' : 49,
            'Amel' : 50,
            'Dmel' : 51,
            'Gacu' : 52,
            'Drer' : 53,
            'Ggal' : 54,
            'Acar' : 55,
            'Oana' : 56,
            'Mdom' : 57,
            'Mmus' : 58,
            'Pabe' : 59,
            'Ptro' : 60,
            'Hsap' : 61}
    return dico

def createListSp():
    Liste = ['Aaeo','Mtub','Tthe','Caur','Mpne','Saur','Bsub','Efae','Spne','Ctra',
            'Bbur','Cjej','Mxan','Gsul','Babo','Wend','Apha','Nmen','Lpne','Ftul',
            'Vcho','Hinf','Ypes','Ecol','Csym','Ckor','Nequ','Paer','Ssol','Taci',
            'Phor','Msmi','Aful','Hbut','Hsal','Mace','Lmaj','Crei','Ppat','Osat',
            'Slyc','Vvin','Atha','Ddis','Spom','Scer','Anid','Ncra','Cele','Amel',
            'Dmel','Gacu','Drer','Ggal','Acar','Oana','Mdom','Mmus','Pabe','Hsap','Ptro']
    return Liste
