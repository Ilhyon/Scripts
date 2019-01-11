#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import os
from pprint import pprint

# ~ listSp = ("apis_mellifera","octopus_bimaculoides","pediculus_humanus","stegodyphus_mimosarum","amphimedon_queenslandica","brugia_malayi","mnemiopsis_leidyi","trichoplax_adhaerens","vitis_vinifera","bigelowiella_natans","emiliania_huxleyi","tetrahymena_thermophila","campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819","brucella_abortus_bv_1_str_9_941","pseudomonas_aeruginosa_mpao1_p2","yersinia_pestis_biovar_microtus_str_91001","anaplasma_phagocytophilum_str_hz","chloroflexus_aurantiacus_j_10_fl","haemophilus_influenzae_rd_kw20","legionella_pneumophila_str_paris","vibrio_cholerae_o1_biovar_el_tor_str_n16961","staphylococcus_aureus_subsp_aureus_n315","francisella_tularensis_subsp_tularensis_schu_s4","escherichia_coli_str_k_12_substr_mg1655","bacillus_subtilis_subsp_subtilis_str_168","mycobacterium_tuberculosis_h37rv","enterococcus_faecalis_v583","mycoplasma_pneumoniae_m129","borrelia_burgdorferi_B31","thermus_thermophilus_HB8","geobacter_sulfurreducens_PCA","wolbachia_endosymbiont_of_drosophila_melanogaster","aquifex_aeolicus_VF5","gardnerella_vaginalis_0288E","myxococcus_xanthus_DK_1622","neisseria_meningitidis_Z2491","pyrobaculum_aerophilum_str_im2","methanobrevibacter_smithii_atcc_35061","cenarchaeum_symbiosum_a","thermoplasma_acidophilum_DSM_1728","methanosarcina_acetivorans_C2A","methanococcus_maripaludis_S2","aeropyrum_pernix_K1","archaeoglobus_fulgidus_DSM_4304","candidatus_korarchaeum_cryptofilum_OPF8")
listSp = ("apis_mellifera", "apis_mellifera")
for sp in listSp :
	print sp
	words = sp.split("_")
	letters = [word[0] for word in words]
	ini = "".join(letters)
	ini = ini.upper()
	filename = "/home/anais/Documents/Data/Genomes/"+sp+"/"+ini+"_transcript_unspliced.txt"
	exists = os.path.isfile(filename)
	if exists :	
		dicoInfo = {}
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					words = l.split('\t')
					if len(words) != 12:
						print words
					idTr = words[1]
					idGene = words[0]
					chrm = words[2]
					Bt = words[3]
					if words[4] or words[5]:
						fUTRs = words[4]
						fUTRe = words[5]
					else :
						fUTRs = ""
						fUTRe = ""
					if words[6] or words[7]:
						tUTRs = words[6]
						tUTRe = words[7]
					else :
						tUTRs = ""
						tUTRe = ""
					if words[8] or words[9] or words[10] :
						exonS = words[8]
						exonE = words[9]
						rank = words[10]
					else :
						exonS = ""
						exonE = ""
						rank = ""
					strand = words[11]
					if (fUTRs or fUTRe) and (exonS or exonE or rank) :
						if idTr not in dicoInfo:
							exonName = "Exon " + rank
							dicoInfo.update({idTr : {"cpt" : 1, exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							dicoInfo.update({idTr : {"cpt" : 0, "5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
						else :
							exonName = "Exon " + rank
							dicoInfo[idTr].update({exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							dicoInfo[idTr]["cpt"] += 1
							dicoInfo[idTr].update({"5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
					elif (tUTRs or tUTRe) and (exonS or exonE or rank):
						if idTr not in dicoInfo :
							exonName = "Exon " + rank
							dicoInfo.update({idTr : {"cpt" : 1, exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							dicoInfo[idTr].update({"3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
						else :
							exonName = "Exon " + rank
							dicoInfo[idTr].update({exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							dicoInfo[idTr]["cpt"] += 1
							dicoInfo[idTr].update({"3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
					else :
						if idTr not in dicoInfo :
							if exonS or exonE or rank :
								exonName = "Exon " + rank
								dicoInfo.update({idTr : {"cpt" : 1, exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							elif tUTRs or tUTRe :
								dicoInfo.update({idTr : {"cpt" : 0, "3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							elif fUTRs or fUTRe :
								dicoInfo.update({idTr : {"cpt" : 0, "5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}}})
							else :
								"Ah"
						else :
							if exonS or exonE or rank :
								exonName = "Exon " + rank
								dicoInfo[idTr].update({exonName : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
								dicoInfo[idTr]["cpt"] += 1
							elif tUTRs or tUTRe :
								dicoInfo[idTr].update({"3UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							elif fUTRs or fUTRe :
								dicoInfo[idTr].update({"5UTR" : {"Gene" : idGene, "Chromosome" : chrm, "Biotype" : Bt, "5UTRstart" : fUTRs, "5UTRend" : fUTRe, "3UTRstart" : tUTRs, "3UTRend" : tUTRe, "exonS" : exonS, "exonE" : exonE, "Rank" : rank, "Strand" : strand}})
							else :
								"Ah"
		words = sp.split("_")
		letters = [word[0] for word in words]
		ini = "".join(letters)
		ini = ini.upper()
		output = open("/home/anais/Documents/Data/Genomes/"+sp+"/"+ini+"_transcript_unspliced.txt","w") # file opening for reading
		for tr in dicoInfo :
			for feature in dicoInfo[tr]:
				if feature != "cpt" and feature != "3UTR" and feature != "5UTR" :
					if "3UTR" in dicoInfo[tr] and "5UTR" in dicoInfo[tr] :
						if feature == "Exon "+str(dicoInfo[tr]["cpt"]) :
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t"+str(dicoInfo[tr]["3UTR"]["3UTRstart"])+"\t"+str(dicoInfo[tr]["3UTR"]["3UTRend"])+"\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
						if feature == "Exon 1" :
							output.write(dicoInfo[tr]["Exon 1"]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr]["Exon 1"]["Chromosome"]+"\t"+dicoInfo[tr]["Exon 1"]["Biotype"]+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRstart"])+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRend"])+"\t\t\t"+str(dicoInfo[tr]["Exon 1"]["exonS"])+"\t"+str(dicoInfo[tr]["Exon 1"]["exonE"])+"\t"+str(dicoInfo[tr]["Exon 1"]["Rank"])+"\t"+dicoInfo[tr]["Exon 1"]["Strand"]+"\n")
						if feature !=  "Exon 1" and feature != "Exon "+str(dicoInfo[tr]["cpt"]):
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t\t\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
					elif "3UTR" in dicoInfo[tr] :
						if feature == "Exon "+str(dicoInfo[tr]["cpt"]) :
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t"+str(dicoInfo[tr]["3UTR"]["3UTRstart"])+"\t"+str(dicoInfo[tr]["3UTR"]["3UTRend"])+"\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
						if feature != "Exon "+str(dicoInfo[tr]["cpt"]):
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t\t\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
					elif "5UTR" in dicoInfo[tr] :
						if feature == "Exon 1" :
							output.write(dicoInfo[tr]["Exon 1"]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr]["Exon 1"]["Chromosome"]+"\t"+dicoInfo[tr]["Exon 1"]["Biotype"]+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRstart"])+"\t"+str(dicoInfo[tr]["5UTR"]["5UTRend"])+"\t\t\t"+str(dicoInfo[tr]["Exon 1"]["exonS"])+"\t"+str(dicoInfo[tr]["Exon 1"]["exonE"])+"\t"+str(dicoInfo[tr]["Exon 1"]["Rank"])+"\t"+dicoInfo[tr]["Exon 1"]["Strand"]+"\n")
						if feature !=  "Exon 1":
							output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+dicoInfo[tr][feature]["Chromosome"]+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t\t\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+dicoInfo[tr][feature]["Strand"]+"\n")
					else:
						output.write(dicoInfo[tr][feature]["Gene"]+"\t"+tr+"\t"+str(dicoInfo[tr][feature]["Chromosome"])+"\t"+dicoInfo[tr][feature]["Biotype"]+"\t\t\t\t\t"+str(dicoInfo[tr][feature]["exonS"])+"\t"+str(dicoInfo[tr][feature]["exonE"])+"\t"+str(dicoInfo[tr][feature]["Rank"])+"\t"+str(dicoInfo[tr][feature]["Strand"])+"\n")
		output.close()
# ~ pprint(dicoInfo)
















