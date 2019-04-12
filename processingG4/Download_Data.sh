#!/bin/bash

# creating all file
for i in anaplasma_phagocytophilum_str_hz anolis_carolinensis apis_mellifera aquifex_aeolicus_vf5 arabidopsis_thaliana archaeoglobus_fulgidus_dsm_4304 aspergillus_nidulans bacillus_subtilis_subsp_subtilis_str_168 borrelia_burgdorferi_b31 brucella_abortus_bv_1_str_9_941 caenorhabditis_elegans campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819 candidatus_korarchaeum_cryptofilum_opf8 cenarchaeum_symbiosum_a chlamydia_trachomatis_d_uw_3_cx chlamydomonas_reinhardtii chloroflexus_aurantiacus_j_10_fl chondrus_crispus danio_rerio dictyostelium_discoideum drosophila_melanogaster emiliania_huxleyi enterococcus_faecalis_v583 escherichia_coli_str_k_12_substr_mg1655 francisella_tularensis_subsp_tularensis_schu_s4 ftp.ensemblgenomes.org gallus_gallus gasterosteus_aculeatus geobacter_sulfurreducens_pca haemophilus_influenzae_rd_kw20 halobacterium_salinarum_r1 homo_sapiens hyperthermus_butylicus_dsm_5456 legionella_pneumophila_str_paris leishmania_major methanobrevibacter_smithii_atcc_35061 methanosarcina_acetivorans_c2a monodelphis_domestica mus_musculus mycobacterium_tuberculosis_h37rv mycoplasma_pneumoniae_m129 myxococcus_xanthus_dk_1622 nanoarchaeum_equitans_kin4_m neisseria_meningitidis_z2491 neurospora_crassa ornithorhynchus_anatinus oryza_sativa pan_troglodyte physcomitrella_patens pongo_abelii pyrobaculum_aerophilum_str_im2 pyrococcus_horikoshii_ot3 saccharomyces_cerevisiae schizosaccharomyces_pombe solanum_lycopersicum staphylococcus_aureus_subsp_aureus_n315 streptococcus_pneumoniae_tigr4 sulfolobus_solfataricus_p2 thermoplasma_acidophilum_dsm_1728 thermus_thermophilus_hb8 vibrio_cholerae_o1_biovar_el_tor_str_n16961 vitis_vinifera wolbachia_endosymbiont_of_drosophila_melanogaster yersinia_pestis_biovar_microtus_str_91001; do mkdir --parents ~/Documents/Data/Genomes/$i;done
for i in $(ls);do mkdir --parents ~/Documents/Data/Genomes/$i/Fasta;done

# retrive sequences
# same command but with different urt because of the different db
# there are some that are from ensembl archive because of the assembly and of the patch assembly
# need to delete file we don't want after
#ensembl
for i in pongo_abelii monodelphis_domestica anolis_carolinensis ornithorhynchus_anatinus danio_rerio gasterosteus_aculeatus;do wget -P ~/Documents/Data/Genomes/$i/Fasta ftp://ftp.ensembl.org/pub/release-95/fasta/$i/dna/*dna.chromosome*.fa.gz;done
wget -P ~/Documents/Data/Genomes/gasterosteus_aculeatus/Fasta ftp://ftp.ensembl.org/pub/release-95/fasta/gasterosteus_aculeatus/dna/*dna.group*.fa.gz;
wget -P ~/Documents/Data/Genomes/homo_sapiens/Fasta ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/*dna.chromosome*.fa.gz
wget -P ~/Documents/Data/Genomes/mus_musculus/Fasta ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/*dna.chromosome*.fa.gz
wget -P ~/Documents/Data/Genomes/gallus_gallus/Fasta  ftp://ftp.ensembl.org/pub/release-75/fasta/gallus_gallus/dna/*dna.chromosome*.fa.gz
wget -P ~/Documents/Data/Genomes/pan_troglodytes/Fasta ftp://ftp.ensembl.org/pub/release-75/fasta/pan_troglodytes/dna/*dna.chromosome*.fa.gz
#metazoa
for i in caenorhabditis_elegans drosophila_melanogaster apis_mellifera; do wget -P ~/Documents/Data/Genomes/$i/Fasta ftp://ftp.ensemblgenomes.org/pub/release-42/metazoa/fasta/$i/dna/*dna.chromosome*.fa.gz;done
#plants
for i in oryza_sativa chlamydomonas_reinhardtii chondrus_crispus physcomitrella_patens solanum_lycopersicum vitis_vinifera arabidopsis_thaliana; do wget -P ~/Documents/Data/Genomes/$i/Fasta ftp://ftp.ensemblgenomes.org/pub/release-42/plants/fasta/$i/dna/*dna.chromosome*.fa.gz;done
#fungi
for i in aspergillus_nidulans neurospora_crassa saccharomyces_cerevisiae schizosaccharomyces_pombe; do wget -P ~/Documents/Data/Genomes/$i/Fasta ftp://ftp.ensemblgenomes.org/pub/release-42/fungi/fasta/$i/dna/*dna.chromosome*.fa.gz;done
#protist
for i in dictyostelium_discoideum emiliania_huxleyi leishmania_major;do wget -P ~/Documents/Data/Genomes/$i/Fasta ftp://ftp.ensemblgenomes.org/pub/release-42/protists/fasta/$i/dna/*dna.chromosome*.fa.gz; done
#bacteria
for i in $(ls);do wget -P ~/Documents/Data/Genomes/$i/Fasta ftp://ftp.ensemblgenomes.org/pub/release-42/bacteria//fasta/bacteria_0_collection/$i/dna/*dna.chromosome*.fa.gz;done

# retrive gtf file
#ensembl
for i in pongo_abelii monodelphis_domestica anolis_carolinensis ornithorhynchus_anatinus danio_rerio gasterosteus_aculeatus; do wget -P ~/Documents/Data/Genomes/$i/ ftp://ftp.ensembl.org/pub/release-95/gtf/$i/*.chr.gtf.gz;done
wget -P ~/Documents/Data/Genomes/homo_sapiens/ ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/*.gtf.gz
wget -P ~/Documents/Data/Genomes/mus_musculus/ ftp://ftp.ensembl.org/pub/release-72/gtf/mus_musculus/*.gtf.gz
wget -P ~/Documents/Data/Genomes/gallus_gallus/ ftp://ftp.ensembl.org/pub/release-75/gtf/gallus_gallus/*.gtf.gz
wget -P ~/Documents/Data/Genomes/pan_troglodytes/ ftp://ftp.ensembl.org/pub/release-75/gtf/pan_troglodytes/*.gtf.gz
#metazoa
for i in caenorhabditis_elegans drosophila_melanogaster apis_mellifera; do wget -P ~/Documents/Data/Genomes/$i/ ftp://ftp.ensemblgenomes.org/pub/release-42/metazoa/gtf/$i/*gtf.gz;done
#plants
for i in oryza_sativa chlamydomonas_reinhardtii chondrus_crispus physcomitrella_patens solanum_lycopersicum vitis_vinifera arabidopsis_thaliana; do wget -P ~/Documents/Data/Genomes/$i ftp://ftp.ensemblgenomes.org/pub/release-42/plants/gtf/$i/*gtf.gz;done
#fungi
for i in aspergillus_nidulans neurospora_crassa saccharomyces_cerevisiae schizosaccharomyces_pombe; do wget -P ~/Documents/Data/Genomes/$i/ ftp://ftp.ensemblgenomes.org/pub/release-42/fungi/gtf/$i/*.gtf.gz;done
#protist
for i in dictyostelium_discoideum emiliania_huxleyi leishmania_major;do wget -P ~/Documents/Data/Genomes/$i/ ftp://ftp.ensemblgenomes.org/pub/release-42/protists/gtf/$i/*.gtf.gz; done
#bacteria
for i in $(ls);do wget -P ~/Documents/Data/Genomes/$i/ ftp://ftp.ensemblgenomes.org/pub/release-42/bacteria//gtf/bacteria_0_collection/$i/*.gtf.gz;done

# miss extract
# miss rename
# miss rm .gz
