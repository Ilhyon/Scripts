#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

def dicoSp():
    {'Homo sapiens' : 'Eukaryota',
    'Pan troglodytes' : 'Eukaryota',
    'Pongo abelii' : 'Eukaryota',
    'Mus musculus' : 'Eukaryota',
    'Monodelphis domestica' : 'Eukaryota',
    'Anolis carolinensis' : 'Eukaryota',
    'Ornithorhynchus anatinus' : 'Eukaryota',
    'Gallus gallus' : 'Eukaryota',
    'Danio rerio' : 'Eukaryota',
    'Gasterosteus aculeatus' : 'Eukaryota',
    'Caenorhabditis elegans' : 'Eukaryota',
    'Drosophila melanogaster' : 'Eukaryota',
    'Apis mellifera' : 'Eukaryota',
    'Oryza sativa' : 'Eukaryota',
    'Chlamydomonas reinhardtii' : 'Eukaryota',
    'Physcomitrella patens' : 'Eukaryota',
    'Solanum lycopersicum' : 'Eukaryota',
    'Vitis vinifera' : 'Eukaryota',
    'Arabidopsis thaliana' : 'Eukaryota',
    'Aspergillus nidulans' : 'Eukaryota',
    'Neurospora crassa' : 'Eukaryota',
    'Saccharomyces cerevisiae' : 'Eukaryota',
    'Schizosaccharomyces pombe' : 'Eukaryota',
    'Dictyostelium discoideum' : 'Eukaryota',
    'Leishmania major' : 'Eukaryota',
    'Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819 strain NCTC 11168' : 'Bacteria',
    'Brucella abortus bv. 1 str. 9-941' : 'Bacteria',
    'Yersinia pestis biovar microtus str 91001' : 'Bacteria',
    'Anaplasma phagocytophilum str. HZ' : 'Bacteria',
    'Chloroflexus aurantiacus J-10-fl' : 'Bacteria',
    'Haemophilus influenzae Rd KW20' : 'Bacteria',
    'Legionella pneumophila str. Paris' : 'Bacteria',
    'Vibrio cholerae O1 biovar El Tor str. N16961' : 'Bacteria',
    'Staphylococcus aureus subsp. aureus N315' : 'Bacteria',
    'Francisella tularensis subsp. tularensis SCHU S4' : 'Bacteria',
    'Escherichia coli str. K-12 substr. MG1655' : 'Bacteria',
    'Bacillus subtilis subsp. subtilis str. 168' : 'Bacteria',
    'Mycobacterium tuberculosis H37Rv' : 'Bacteria',
    'Enterococcus faecalis V583' : 'Bacteria',
    'Mycoplasma pneumoniae M129' : 'Bacteria',
    'Streptococcus pneumoniae TIGR4' : 'Bacteria',
    'Borreliella burgdorferi B31' : 'Bacteria',
    'Thermus thermophilus HB8' : 'Bacteria',
    'Geobacter sulfurreducens PCA' : 'Bacteria',
    'Wolbachia endosymbiont of drosophila melanogaster' : 'Bacteria',
    'Aquifex aeolicus VF5' : 'Bacteria',
    'Myxococcus xanthus DK 1622' : 'Bacteria',
    'Neisseria meningitidis Z2491' : 'Bacteria',
    'Chlamydia trachomatis D/UW-3/CX' : 'Bacteria',
    'Nanoarchaeum equitans Kin4-M' : 'Archaea',
    'Pyrobaculum aerophilum str. IM2' : 'Archaea',
    'Methanobrevibacter smithii ATCC 35061' : 'Archaea',
    'Cenarchaeum symbiosum A' : 'Archaea',
    'Saccharolobus solfataricus P2' : 'Archaea',
    'Thermoplasma acidophilum DSM 1728' : 'Archaea',
    'Methanosarcina acetivorans C2A' : 'Archaea',
    'Pyrococcus horikoshii OT3' : 'Archaea',
    'Archaeoglobus fulgidus dsm 4304' : 'Archaea',
    'Candidatus Korarchaeum cryptofilum OPF8' : 'Archaea',
    'Halobacterium salinarum R1' : 'Archaea',
    'Hyperthermus butylicus DSM 5456' : 'Archaea'}

def getDicoSpId():
    dico = {'APSH' : 'anaplasma_phagocytophilum_str_hz',
            'AC' : 'anolis_carolinensis',
            'AM' : 'apis_mellifera',
            'AAV' : 'aquifex_aeolicus_vf5',
            'AT' : 'arabidopsis_thaliana',
            'AFD4' : 'archaeoglobus_fulgidus_dsm_4304',
            'AN' : 'aspergillus_nidulans',
            'BSSSS1' : 'bacillus_subtilis_subsp_subtilis_str_168',
            'BBB' : 'borrelia_burgdorferi_b31',
            'BAB1S99' : 'brucella_abortus_bv_1_str_9_941',
            'CE' : 'caenorhabditis_elegans',
            'CJSJN1A7' : 'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819',
            'CKCO' : 'candidatus_korarchaeum_cryptofilum_opf8',
            'CSA' : 'cenarchaeum_symbiosum_a',
            'CTDU3C' : 'chlamydia_trachomatis_d_uw_3_cx',
            'CR' : 'chlamydomonas_reinhardtii',
            'CAJ1F' : 'chloroflexus_aurantiacus_j_10_fl',
            'CC' : 'chondrus_crispus',
            'DR' : 'danio_rerio',
            'DD' : 'dictyostelium_discoideum',
            'DM' : 'drosophila_melanogaster',
            'EFV' : 'enterococcus_faecalis_v583',
            'ECSK1SM' : 'escherichia_coli_str_k_12_substr_mg1655',
            'FTSTSS' : 'francisella_tularensis_subsp_tularensis_schu_s4',
            'GG' : 'gallus_gallus',
            'GA' : 'gasterosteus_aculeatus',
            'GSP' : 'geobacter_sulfurreducens_pca',
            'HIRK' : 'haemophilus_influenzae_rd_kw20',
            'HSR' : 'halobacterium_salinarum_r1',
            'HS' : 'homo_sapiens',
            'HBD5' : 'hyperthermus_butylicus_dsm_5456',
            'LPSP' : 'legionella_pneumophila_str_paris',
            'LM' : 'leishmania_major',
            'MSA3' : 'methanobrevibacter_smithii_atcc_35061',
            'MAC' : 'methanosarcina_acetivorans_c2a',
            'MD' : 'monodelphis_domestica',
            'MM' : 'mus_musculus',
            'MTH' : 'mycobacterium_tuberculosis_h37rv',
            'MPM' : 'mycoplasma_pneumoniae_m129',
            'MXD1' : 'myxococcus_xanthus_dk_1622',
            'NEKM' : 'nanoarchaeum_equitans_kin4_m',
            'NMZ' : 'neisseria_meningitidis_z2491',
            'NC' : 'neurospora_crassa',
            'OA' : 'ornithorhynchus_anatinus',
            'OS' : 'oryza_sativa',
            'PT' : 'pan_troglodytes',
            'PP' : 'physcomitrella_patens',
            'PA' : 'pongo_abelii',
            'PASI' : 'pyrobaculum_aerophilum_str_im2',
            'PHO' : 'pyrococcus_horikoshii_ot3',
            'SC' : 'saccharomyces_cerevisiae',
            'SP' : 'schizosaccharomyces_pombe',
            'SL' : 'solanum_lycopersicum',
            'SASAN' : 'staphylococcus_aureus_subsp_aureus_n315',
            'SPT' : 'streptococcus_pneumoniae_tigr4',
            'SSP' : 'sulfolobus_solfataricus_p2',
            'TAD1' : 'thermoplasma_acidophilum_dsm_1728',
            'TTH' : 'thermus_thermophilus_hb8',
            'VCOBETSN' : 'vibrio_cholerae_o1_biovar_el_tor_str_n16961',
            'VV' : 'vitis_vinifera',
            'WEODM' : 'wolbachia_endosymbiont_of_drosophila_melanogaster',
            'YPBMS9' : 'yersinia_pestis_biovar_microtus_str_91001'}
    return dico