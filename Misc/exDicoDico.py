#!/usr/bin/python
# -*- coding: utf-8 -*-

####################################################
# Lecture du jeu d'apprentissage pour construire   #
# la table de probabilité                          #
####################################################

import argparse
import re
from pprint import pprint
 
# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='k-nn classifier.')
parser.add_argument('--training', required=True, help='File containing training examples (class,att1,...')
parser.add_argument('--sample', required=False, help='File containing new objects to be classified.')
parser.add_argument('--k', required=False, default=5, help='Number of nearest neighbor to consider (default 5)')
parser.add_argument('--normalize', required=False, action="store_true", help='Normalize values.')
opt = parser.parse_args()

# AFFICHAGE DES PARAMETRES PASSES EN LIGNE DE COMMANDEcla	
# print opt

import csv
 
# LOAD TRAINING DATA
y = []
data = []
occurences = {}
att = []
with open(opt.training) as csvfile: # ouvre le fichier spécifié après --training
	csvreader = csv.reader(csvfile, delimiter=';') # on précise que le séparateur est une virgule
	head = ("Query_ID","Subject_ID","Domain_Structure_Subject","FD_Name_Subject","Domain_Structure_Query")
	cpt_DS = 0 # pour connaitre le nombre tot d'élément dans le jeu de donnée training
	cpt_FDN = 0
	for row in csvreader: # read each line
		c = row[0] # Query id
		a = head[2] # domain structure
		b = head[3] # "Fonctionnal domain"
		v = row[2] # domain structure (MSD, SBP ...)
		x = row[3] # FD_Name (pfam, cog...)
		z = row[4] # domaine vrai
		if c not in occurences: #si le gène n'existe pas dans le dico, on le crée
			occurences[c] = {}
			occurences[c]["Domain_vrai"] = z
		if a not in occurences[c]: # si l'attribut n'exsite pas
			occurences[c][a] = {} # on le crée (dico avec les modalités)
		if v not in occurences[c][a]:# si domain structure (MSD, SBP ...) n'existe pas
			occurences[c][a][v] = {} # on la crée et elle commence à 1
			occurences[c][a][v]["cpt_DS"] = 0 # cpt qui compte le nombre de DS différents
		if b not in occurences[c][a][v]:
			occurences[c][a][v][b] = {} # on le crée (dico avec les modalités)
		if x not in occurences[c][a][v][b]:
			occurences[c][a][v][b][x] = 1
			occurences[c][a][v]["cpt_DS"] += 1
			#occurences[c]["cpt"] += 1
		else: # la modalité existe donc on ajoute 1
			#occurences[c]["cpt"] += 1
			occurences[c][a][v][b][x] += 1
			occurences[c][a][v]["cpt_DS"] += 1


cpt = 0
tmpCpt = 0
qi = "truc"
tmpDom = "truc"
print("Query_ID ; Domain_Predit ; Pourcentage ; Domain_vrai_QI")
for c in occurences:
	for a in occurences[c]:
		if a == "Domain_vrai":
			dv = occurences[c][a]
			qi = c
		if a != "Domain_vrai" and qi == c:
			for v in occurences[c][a]:
				if occurences[c][a][v]["cpt_DS"] > tmpCpt :
					tmpCpt = occurences[c][a][v]["cpt_DS"]
					tmpDom = v
				cpt += occurences[c][a][v]["cpt_DS"]
			name = "% max = " + tmpDom
			print(c, ";", v, ";", (tmpCpt / cpt) * 100, ";", dv, "")
			cpt = 0
			tmpCpt = 0

#print(cpt)
