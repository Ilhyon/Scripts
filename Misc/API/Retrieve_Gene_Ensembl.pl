use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-eg-publicsql.ebi.ac.uk',
    -port => 4157,
);

my @list_species = ("Homo sapiens","Pan troglodytes","Pongo abelii","Mus musculus","Monodelphis domestica","Anolis carolinensis","Ornithorhynchus anatinus","Gallus gallus","Danio rerio","Gasterosteus aculeatus","Xenopus tropicalis","Ciona savignyi","Caenorhabditis elegans","Drosophila melanogaster","Anopheles gambiae","Apis mellifera","Octopus bimaculoides","Pediculus humanus","Stegodyphus mimosarum","Amphimedon queenslandica","Brugia malayi","Mnemiopsis leidyi","Trichoplax adhaerens","Chlamydomonas reinhardtii","Chondrus crispus","Amborella trichopoda","Oryza sativa Japonica Group","Physcomitrella patens","Selaginella moellendorffii","Solanum lycopersicum","Vitis vinifera","Arabidopsis thaliana","Aspergillus nidulans","Neurospora crassa","Saccharomyces Cerevisiae","Schizosaccharomyces pombe","Bigelowiella natans","Dictyostelium discoideum (om)","Emiliania huxleyi (très nombreux sur terre)","Leishmania major (maladie)","Tetrahymena thermophila","Thermoanaerobacter kivui","Phytoplasma onion yellows","Mycoplasma mycoides","Mycoplasma pneumoniae","Staphylococcus aureus","Listeria monocytogenes ","Bacillus subtilis","Enterococcus faecalis","Streptococcus pneumoniae","Chlamydia trachomatis","Borrelia burgdorferi","Mycobacterium leprae","Mycobacterium tuberculosis","Thermus thermophilus","Synechococcus sp. WH8102","Geobacter sulfurreducens","Campylobacter jejuni","Wolbachia","Brucella abortus","Nitrosomonas communis","Pseudomonas aeruginosa","Pseudomonas putida","Escherichia coli","Yersinia pestis","Anaplasma phagocytophilum","Aquifex aeolicus","Caulobacter crescentus","Cenarchaeum symbiosum","Chloroflexus aurantiacus","Coxiella burnetii","Enterobacter cloacae","Francisella tularensis","Gardnerella vaginalis","Haemophilus influenzae","Helicobacter pylori","Klebsiella pneumoniae","Lactobacillus plantarum","Legionella pneumophila","Leuconostoc mesenteroides","Moraxella catarrhalis","Myxococcus xanthus","Neisseria meningitidis","Vibrio cholerae","Nanoarchaeum equitans","Pyrobaculum aerophilum","Sulfolobus solfataricus","Thermoplasma acidophilum","Methanosarcina acetivorans","Pyrococcus horikoshii","Methanopyrus kandleri","Methanococcus maripaludis","Aeropyrum pernix","Archaeoglobus fulgidus","Candidatus Korarchaeum","Haloarcula marismortui ","Halobacterium salinarum","Haloferax volcanii","Hyperthermus butylicus","Methanobrevibacter","Nanoarchaeum equitans");

my $pan_homology_adaptor  = $registry->get_adaptor('pan_homology', 'compara', 'GeneMember');
my @gene_members = $pan_homology_adaptor->fetch_all();
#~ print @gene_members; ARRAY(0x1ac7b80)
#~ my $test = $pan_homology_adaptor->fetch_all();
#~ print $test; ARRAY(0xbeb7bf88)

foreach my $gene (@gene_members) {
  #~ print $gene; # array
  foreach my $i (@{$gene}){
	  my $genome_DB_id = $i->genome_db_id();
	  print $genome_DB_id, "\n";
	  }
}

