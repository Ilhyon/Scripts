#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import argparse
import networkx as nx

def iniListSpecie():
    list = ['homo_sapiens','pan_troglodytes','pongo_abelii','mus_musculus',
            'monodelphis_domestica','anolis_carolinensis',
            'ornithorhynchus_anatinus','gallus_gallus','danio_rerio',
            'gasterosteus_aculeatus','caenorhabditis_elegans',
            'drosophila_melanogaster','apis_mellifera','oryza_sativa',
            'chlamydomonas_reinhardtii','physcomitrella_patens',
            'solanum_lycopersicum','vitis_vinifera','arabidopsis_thaliana',
            'aspergillus_nidulans','neurospora_crassa','saccharomyces_cerevisiae',
            'schizosaccharomyces_pombe','dictyostelium_discoideum',
            'emiliania_huxleyi','leishmania_major',
            'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819',
            'brucella_abortus_bv_1_str_9_941','yersinia_pestis_biovar_microtus_str_91001',
            'anaplasma_phagocytophilum_str_hz','chloroflexus_aurantiacus_j_10_fl',
            'haemophilus_influenzae_rd_kw20','legionella_pneumophila_str_paris',
            'vibrio_cholerae_o1_biovar_el_tor_str_n16961',
            'staphylococcus_aureus_subsp_aureus_n315',
            'francisella_tularensis_subsp_tularensis_schu_s4',
            'escherichia_coli_str_k_12_substr_mg1655','bacillus_subtilis_subsp_subtilis_str_168',
            'mycobacterium_tuberculosis_h37rv','enterococcus_faecalis_v583',
            'mycoplasma_pneumoniae_m129','streptococcus_pneumoniae_tigr4',
            'borrelia_burgdorferi_b31','thermus_thermophilus_hb8','geobacter_sulfurreducens_pca',
            'wolbachia_endosymbiont_of_drosophila_melanogaster','aquifex_aeolicus_vf5',
            'myxococcus_xanthus_dk_1622','neisseria_meningitidis_z2491',
            'chlamydia_trachomatis_d_uw_3_cx','nanoarchaeum_equitans_kin4_m',
            'pyrobaculum_aerophilum_str_im2','methanobrevibacter_smithii_atcc_35061',
            'cenarchaeum_symbiosum_a','sulfolobus_solfataricus_p2',
            'thermoplasma_acidophilum_dsm_1728','methanosarcina_acetivorans_c2a',
            'pyrococcus_horikoshii_ot3','archaeoglobus_fulgidus_dsm_4304',
            'candidatus_korarchaeum_cryptofilum_opf8','halobacterium_salinarum_r1',
            'hyperthermus_butylicus_dsm_5456']
    return list

def importHomology(filename):
    G = nx.Graph()
    spList = iniListSpecie()
    with open(filename) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines: # browse all lines
            if not l.startswith('gene_stable_id') and l:
                words = l.split('\t')
                gene1 = words[0]
                gene2 = words[5]
                homologyT = words[4]
                specie1 = words[2]
                specie2 = words[7]
                if specie1 in spList and specie2 in specie2:
                    G.add_node(gene1, specie=specie1)
                    G.add_node(gene2, specie=specie2)
                    G.add_edge(gene1, gene2, homology=homologyT)

def main(filename):
    importHomology(filename)

if __name__ == '__main__':
    filename = "/home/anais/Documents/Data/Homology/Compara_default_homologies.tsv"
    main(filename)
