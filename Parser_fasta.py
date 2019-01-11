#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import os
import math
from pprint import pprint

# ~ listSp = ("pyrobaculum_aerophilum_str_im2","methanobrevibacter_smithii_atcc_35061","staphylococcus_aureus_subsp_aureus_n315","francisella_tularensis_subsp_tularensis_schu_s4","escherichia_coli_str_k_12_substr_mg1655","bacillus_subtilis_subsp_subtilis_str_168","mycobacterium_tuberculosis_h37rv","enterococcus_faecalis_v583","mycoplasma_pneumoniae_m129","streptococcus_pneumoniae_TIGR4","borrelia_burgdorferi_B31","thermus_thermophilus_HB8","geobacter_sulfurreducens_PCA","wolbachia_endosymbiont_of_drosophila_melanogaster","aquifex_aeolicus_VF5","cenarchaeum_symbiosum_a","gardnerella_vaginalis_0288E","myxococcus_xanthus_DK_1622","neisseria_meningitidis_Z2491","sulfolobus_solfataricus_P2","thermoplasma_acidophilum_DSM_1728","methanosarcina_acetivorans_C2A","pyrococcus_horikoshii_OT3","methanococcus_maripaludis_S2","aeropyrum_pernix_K1","archaeoglobus_fulgidus_DSM_4304","candidatus_korarchaeum_cryptofilum_OPF8","halobacterium_salinarum_R1","campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819","brucella_abortus_bv_1_str_9_941","pseudomonas_aeruginosa_mpao1_p2","yersinia_pestis_biovar_microtus_str_91001","anaplasma_phagocytophilum_str_hz","chloroflexus_aurantiacus_j_10_fl","haemophilus_influenzae_rd_kw20","legionella_pneumophila_str_paris","vibrio_cholerae_o1_biovar_el_tor_str_n16961","chlamydia_trachomatis_d_uw_3_cx","nanoarchaeum_equitans_kin4_m")
listSp = ("apis_mellifera","apis_mellifera")
for sp in listSp :	
	direcoty = "/home/anais/Documents/Data/Genomes/"+sp+"/"
	dicoChromosome = {}
	listFile = os.listdir(direcoty)
	for filename in listFile :
		with open(filename) as f: # file opening
			content = f.read()
			l = content.split('\n')
			if l[0].startswith('>'):
				header = l[0]
				chrm = header.split(' ')[0][1:]
				sequence = "".join(l[1:])
				dicoChromosome.update({chrm : sequence})
	filename = "/home/anais/Documents/Data/Genomes/"+sp+"/"+sp+".gtf"
	exists = os.path.isfile(filename)
	if exists :	
		dicoGene = {}
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					words=l.split('\t')
					idGene = words[8].split(';')[0].split('"')[1]
					feature = words[2]
					if feature == "gene" :
						chrm = words[0]
						startFeature = words[3]
						endFeature = words[4]
						strand = words[6]
						if strand == "+":
							strand = 1
						elif strand == "-":
							strand = -1
						geneSequence = dicoChromosome[chrm][(int(startFeature)-1):(int(endFeature)-1+1)]
						dicoGene.update({idGene : {"Chromosome" : chrm, "Start" : startFeature, "End" : endFeature, "Strand" : str(strand), "Sequence" : geneSequence}})
	
	words = sp.split("_")
	letters = [word[0] for word in words]
	ini = "".join(letters)
	ini = ini.upper()
	output = open("/home/anais/Documents/Data/Genomes/"+sp+"/"+ini+"_gene_unspliced.txt","w") # file opening for reading
	for gene in dicoGene :
		chromosome = dicoGene[gene]["Chromosome"]
		start = dicoGene[gene]["Start"]
		end = dicoGene[gene]["End"]
		Sequence = dicoGene[gene]["Sequence"]
		strand = str(dicoGene[gene]["Strand"])
		output.write(">"+gene+"|"+chromosome+"|"+start+"|"+end+"|"+strand+"\n")
		nbLine = math.ceil(float(int(end)-int(start))/60)
		cpt1 = 0
		cpt2 = 60
		if nbLine < 0 :
			print AH
		if strand == "-1" :
			reverse = ""
			for n in Sequence:
				if n == "A" :
					tmp = "T"
				elif n == "T" :
					tmp = "A"
				elif n == "G" :
					tmp = "C"
				elif n == "C" :
					tmp = "G"
				else :
					tmp = n
				reverse = reverse + tmp
			Sequence = reverse[::-1]
		for i in range(0,int(nbLine)) :
			output.write(Sequence[cpt1:cpt2]+"\n")
			cpt1 += 60
			cpt2 += 60
	
# ~ pprint(dicoGene)
