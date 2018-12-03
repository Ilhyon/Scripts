#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-eg-publicsql.ebi.ac.uk',
    -port => 4157,
);

print "Succeed to connect to the DB";

my @list_species = ("homo_sapiens","pan_troglodytes","pongo_abelii","mus_musculus","monodelphis_domestica","anolis_carolinensis","ornithorhynchus_anatinus","gallus_gallus","danio_rerio","gasterosteus_aculeatus","xenopus_tropicalis","ciona_savignyi","caenorhabditis_elegans","drosophila_melanogaster","anopheles_gambiae","apis_mellifera","octopus_bimaculoides","pediculus_humanus","stegodyphus_mimosarum","amphimedon_queenslandica","brugia_malayi","mnemiopsis_leidyi","trichoplax_adhaerens","chlamydomonas_reinhardtii","chondrus_crispus","amborella_trichopoda","physcomitrella_patens","selaginella_moellendorffii","solanum_lycopersicum","vitis_vinifera","arabidopsis_thaliana","aspergillus_nidulans","neurospora_crassa","saccharomyces_cerevisiae","schizosaccharomyces_pombe","bigelowiella_natans","dictyostelium_discoideum","emiliania_huxleyi","leishmania_major","tetrahymena_thermophila","mycoplasma_pneumoniae_m129","streptococcus_pneumoniae_TIGR4","borrelia_burgdorferi_B31","thermus_thermophilus_HB8","geobacter_sulfurreducens_PCA","wolbachia_endosymbiont_of_drosophila_melanogaster","aquifex_aeolicus_VF5","caulobacter_crescentus_CB15","cenarchaeum_symbiosum_a","coxiella_burnetii_rsa_493","gardnerella_vaginalis_0288E","helicobacter_pylori_26695","lactobacillus_plantarum_WCFS1","moraxella_catarrhalis_7169","myxococcus_xanthus_DK_1622","neisseria_meningitidis_Z2491","sulfolobus_solfataricus_P2","thermoplasma_acidophilum_DSM_1728","methanosarcina_acetivorans_C2A","methanopyrus_kandleri_AV19","methanococcus_maripaludis_S2","aeropyrum_pernix_K1","archaeoglobus_fulgidus_DSM_4304","candidatus_korarchaeum_cryptofilum_OPF8","haloarcula_marismortui_ATCC_43049","halobacterium_salinarum_R1","hyperthermus_butylicus_DSM_5456");
my @list_taxonID = ("39947","158879","169963","224308","272561","192222","262698","1131758","511145","229193","212042","324602","716541","177416","71421","272620","297246","203120","243277","228908","178306","204669","420247","228908");

my $panH_genomeDB_adaptor = $registry->get_adaptor('pan_homology', 'compara', 'GenomeDB');

#~ print $panH_genomeDB_adaptor; #Bio::EnsEMBL::Compara::DBSQL::GenomeDBAdaptor=HASH(0x1790b40)

my @genomes_db = ();

# get all the genome db for the list of species (with names)
foreach my $sp (@list_species) {
	my $genome_db = $panH_genomeDB_adaptor->fetch_by_name_assembly($sp); # Bio::EnsEMBL::Compara::GenomeDB
	if($genome_db){
		push @genomes_db, $genome_db;
	}
}

# get all the genome db for the list of species (with taxon ID)
foreach my $taxonID (@list_taxonID) {
	my $genome_db = $panH_genomeDB_adaptor->fetch_all_by_taxon_id($taxonID);
	if($genome_db){
		push @genomes_db, $genome_db;
	}
}

print "Succeed to retrieve the genomes DB";

#~ print scalar @genomes_db; # 91 same number as species given in entry

my $db_adaptor = $genomes_db[1]->db_adaptor();
print $db_adaptor,"\n"; #
my $gene_adaptor = $db_adaptor->get_GeneAdaptor();
print $gene_adaptor,"\n"; #
