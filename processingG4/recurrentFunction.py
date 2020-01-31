#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

def createDicoFamily():
	"""Creates a dictionnary with all classes and subclasses of transcripts.

	:returns: dicoFam, contains all classes and subclasses.
	:rtype: dictionary
	"""
	dicoFam = {'Coding' : ['IG_C_gene', 'IG_D_gene', 'IG_J_gene',
							'IG_LV_gene', 'IG_M_gene', 'IG_V_gene',
							'IG_Z_gene', 'nonsense_mediated_decay',
							'nontranslating_CDS', 'non_stop_decay',
							'protein_coding', 'TR_C_gene', 'TR_D_gene',
							'TR_gene', 'TR_J_gene', 'TR_V_gene'],
				'Pseudogene' : ['transcribed_unitary_pseudogene',
								'disrupted_domain', 'IG_C_pseudogene',
								'IG_J_pseudogene', 'IG_pseudogene',
								'IG_V_pseudogene', 'processed_pseudogene',
								'pseudogene',
								'transcribed_processed_pseudogene',
								'transcribed_unprocessed_pseudogene',
								'translated_processed_pseudogene',
								'translated_unprocessed_pseudogene',
								'TR_J_pseudogene', 'TR_V_pseudogene',
								'unitary_pseudogene', 'unprocessed_pseudogene',
								'polymorphic_pseudogene'],
				'LongNC' : ['macro_lncRNA', 'bidirectional_promoter_lncRNA',
							'sense_intronic', '3prime_overlapping_ncRNA',
							'ambiguous_orf', 'antisense',
							'lincRNA', 'ncrna_host','non_coding',
							'processed_transcript', 'retained_intron',
							'sense_overlapping'],
				'ShortNC' : ['vaultRNA', 'scaRNA', 'miRNA',
							'miRNA_pseudogene', 'misc_RNA',
							'misc_RNA_pseudogene', 'Mt_rRNA',
							'Mt_tRNA', 'Mt_tRNA_pseudogene',
							'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA',
							'RNase_P_RNA', 'rRNA', 'rRNA_pseudogene',
							'scRNA', 'scRNA_pseudogene', 'snlRNA',
							'snoRNA', 'snoRNA_pseudogene', 'snRNA',
							'snRNA_pseudogene', 'SRP_RNA', 'tmRNA',
							'tRNA', 'tRNA_pseudogene','ribozyme'],
				'Predictif' : ['TEC']}
	return dicoFam

def createDicoType():
	"""Creates a dictionnary with all sublclasses and classes of transcripts.

	:returns: dicoFam, contains all classes and subclasses.
	:rtype: dictionary
	"""
	dicoFam = {'IG_C_gene' : 'Coding',
			'IG_D_gene' : 'Coding',
			'IG_J_gene' : 'Coding',
			'IG_LV_gene' : 'Coding',
			'IG_M_gene' : 'Coding',
			'IG_V_gene' : 'Coding',
			'IG_Z_gene' : 'Coding',
			'nonsense_mediated_decay' : 'Coding',
			'nontranslating_CDS' : 'Coding',
			'non_stop_decay' : 'Coding',
			'protein_coding' : 'Coding',
			'TR_C_gene' : 'Coding',
			'TR_D_gene' : 'Coding',
			'TR_gene' : 'Coding',
			'TR_J_gene' : 'Coding',
			'TR_V_gene' : 'Coding',
			'transcribed_unitary_pseudogene' : 'Pseudogene',
			'disrupted_domain' : 'Pseudogene',
			'IG_C_pseudogene' : 'Pseudogene',
			'IG_J_pseudogene' : 'Pseudogene',
			'IG_pseudogene' : 'Pseudogene',
			'IG_V_pseudogene' : 'Pseudogene',
			'processed_pseudogene' : 'Pseudogene',
			'pseudogene' : 'Pseudogene',
			'transcribed_processed_pseudogene' : 'Pseudogene',
			'transcribed_unprocessed_pseudogene' : 'Pseudogene',
			'translated_processed_pseudogene' : 'Pseudogene',
			'translated_unprocessed_pseudogene' : 'Pseudogene',
			'TR_J_pseudogene' : 'Pseudogene',
			'TR_V_pseudogene' : 'Pseudogene',
			'unitary_pseudogene' : 'Pseudogene',
			'unprocessed_pseudogene' : 'Pseudogene',
			'polymorphic_pseudogene' : 'Pseudogene',
			'macro_lncRNA' : 'LongNC',
			'bidirectional_promoter_lncRNA' : 'LongNC',
			'sense_intronic' : 'LongNC',
			'3prime_overlapping_ncRNA' : 'LongNC',
			'ambiguous_orf' : 'LongNC',
			'antisense' : 'LongNC',
			'lincRNA' : 'LongNC',
			'ncrna_host''non_coding' : 'LongNC',
			'processed_transcript' : 'LongNC',
			'retained_intron' : 'LongNC',
			'sense_overlapping' : 'LongNC',
			'vaultRNA' : 'ShortNC',
			'scaRNA' : 'ShortNC',
			'miRNA' : 'ShortNC',
			'miRNA_pseudogene' : 'ShortNC',
			'misc_RNA' : 'ShortNC',
			'misc_RNA_pseudogene' : 'ShortNC',
			'Mt_rRNA' : 'ShortNC',
			'Mt_tRNA' : 'ShortNC',
			'Mt_tRNA_pseudogene' : 'ShortNC',
			'ncRNA' : 'ShortNC',
			'pre_miRNA' : 'ShortNC',
			'RNase_MRP_RNA' : 'ShortNC',
			'RNase_P_RNA' : 'ShortNC',
			'rRNA' : 'ShortNC',
			'rRNA_pseudogene' : 'ShortNC',
			'scRNA' : 'ShortNC',
			'scRNA_pseudogene' : 'ShortNC',
			'snlRNA' : 'ShortNC',
			'snoRNA' : 'ShortNC',
			'snoRNA_pseudogene' : 'ShortNC',
			'snRNA' : 'ShortNC',
			'snRNA_pseudogene' : 'ShortNC',
			'SRP_RNA' : 'ShortNC',
			'tmRNA' : 'ShortNC',
			'tRNA' : 'ShortNC',
			'tRNA_pseudogene' : 'Pseudogene',
			'ribozyme' : 'ShortNC'}
	return dicoFam

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
