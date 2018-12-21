#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
import os
from pprint import pprint

listSp = ("pyrobaculum_aerophilum_str_im2","methanobrevibacter_smithii_atcc_35061","staphylococcus_aureus_subsp_aureus_n315","francisella_tularensis_subsp_tularensis_schu_s4","escherichia_coli_str_k_12_substr_mg1655","bacillus_subtilis_subsp_subtilis_str_168","mycobacterium_tuberculosis_h37rv","enterococcus_faecalis_v583","mycoplasma_pneumoniae_m129","streptococcus_pneumoniae_TIGR4","borrelia_burgdorferi_B31","thermus_thermophilus_HB8","geobacter_sulfurreducens_PCA","wolbachia_endosymbiont_of_drosophila_melanogaster","aquifex_aeolicus_VF5","cenarchaeum_symbiosum_a","gardnerella_vaginalis_0288E","myxococcus_xanthus_DK_1622","neisseria_meningitidis_Z2491","sulfolobus_solfataricus_P2","thermoplasma_acidophilum_DSM_1728","methanosarcina_acetivorans_C2A","pyrococcus_horikoshii_OT3","methanococcus_maripaludis_S2","aeropyrum_pernix_K1","archaeoglobus_fulgidus_DSM_4304","candidatus_korarchaeum_cryptofilum_OPF8","halobacterium_salinarum_R1","campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819","brucella_abortus_bv_1_str_9_941","pseudomonas_aeruginosa_mpao1_p2","yersinia_pestis_biovar_microtus_str_91001","anaplasma_phagocytophilum_str_hz","chloroflexus_aurantiacus_j_10_fl","haemophilus_influenzae_rd_kw20","legionella_pneumophila_str_paris","vibrio_cholerae_o1_biovar_el_tor_str_n16961","chlamydia_trachomatis_d_uw_3_cx","nanoarchaeum_equitans_kin4_m")
#~ listSp = ("pyrobaculum_aerophilum_str_im2")
for sp in listSp :
	filename = "/home/anais/Documents/Data/"+sp+"/"+sp+".gtf"
	exists = os.path.isfile(filename)
	if exists :	
		dicoFeature = {}
		with open(filename) as f: # file opening
			content = f.read()
			lines = content.split('\n')
			for l in lines: # browse all lines
				if not l.startswith('#') and l:
					#~ print l
					words=l.split('\t')
					attribute = words[8].split(';')
					idGene = attribute[0].split('"')[1]
					feature = words[2]
					chrm = words[0]
					startFeature = words[3]
					endFeature = words[4]
					strand = words[6]
					if idGene not in dicoFeature :
						if feature == "gene" :
							biotype = attribute[2].split('"')[1]
							dicoFeature[idGene] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
						elif feature == "transcript" :
							idTr = attribute[1].split('"')[1]
							biotype = attribute[5].split('"')[1]
							dicoFeature[idGene] = {"Transcript" : idTr}
							dicoFeature[idGene]["Transcript"][idTr] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
						elif feature == "exon" :
							rank = attribute[2].split('"')[1]
							idTr = attribute[1].split('"')[1]
							idExon = attribute[-2].split('"')[1]
							biotype = attribute[5].split('"')[1]
							dicoFeature[idGene] = {"Transcript" : idTr}
							dicoFeature[idGene]["Transcript"][idTr] = {"Exon" : idExon}
							dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
					else:
						if feature == "gene" :
							biotype = attribute[2].split('"')[1]
							pprint(dicoFeature[idGene])
							print "Chromosome : "+chrm+" | Start : "+startFeature+" | End :"+endFeature+" | Biotype : "+biotype+" | Strand : "+strand
						elif feature == "transcript" :
							idTr = attribute[1].split('"')[1]
							biotype = attribute[5].split('"')[1]
							if "Transcript" in dicoFeature[idGene]:
								if idTr not in dicoFeature[idGene]["Transcript"] :
									dicoFeature[idGene]["Transcript"][idTr] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
								else :
									biotype = attribute[2].split('"')[1]
									pprint(dicoFeature[idGene]["Transcript"][idTr])
									print "Chromosome : "+chrm+" | Start : "+startFeature+" | End :"+endFeature+" | Biotype : "+biotype+" | Strand : "+strand
							else : 
								dicoFeature[idGene]["Transcript"] = {idTr:{}}
								dicoFeature[idGene]["Transcript"][idTr] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand}
						elif feature == "exon" :
							rank = attribute[2].split('"')[1]
							idTr = attribute[1].split('"')[1]
							biotype = attribute[5].split('"')[1]
							idExon = attribute[-2].split('"')[1]
							if "Transcript" in dicoFeature[idGene]:
								if idTr not in dicoFeature[idGene]["Transcript"] :
									dicoFeature[idGene]["Transcript"][idTr] = {"Exon" : idExon}
									dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
								else :
									if "Exon" in dicoFeature[idGene]["Transcript"][idTr]:
										if idExon not in dicoFeature[idGene]["Transcript"][idTr]["Exon"]:
											dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
										else:	
											pprint(dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon])
											print "Chromosome : "+chrm+" | Start : "+startFeature+" | End :"+endFeature+" | Biotype : "+biotype+" | Strand : "+strand
									else :
										dicoFeature[idGene]["Transcript"][idTr]["Exon"] ={idExon : {}}
										dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
							else : 
								dicoFeature[idGene].update({"Transcript" : idTr})
								dicoFeature[idGene]["Transcript"][idTr] = {"Exon" : idExon}
								dicoFeature[idGene]["Transcript"][idTr]["Exon"][idExon] = {"Chromosome" : chrm, "Start" : startFeature, "End" :endFeature, "Biotype" : biotype, "Strand" : strand, "Rank" : rank}
		pprint(dicoFeature)
	
	#~ else :
		#~ print filename
