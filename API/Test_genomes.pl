use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::Species;
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-eg-publicsql.ebi.ac.uk',
    -port => 4157,
);



my $genome_db_adaptor = $registry->get_adaptor('metazoa', 'compara', 'GenomeDB');
my $genome_db = $genome_db_adaptor->fetch_by_name_assembly("caenorhabditis_elegans"); # Bio::EnsEMBL::Compara::GenomeDB
print $genome_db,"\n";
my $db_adaptor = $genome_db->db_adaptor(); # Bio::EnsEMBL::DBSQL::DBAdaptor
my $gene_adaptor = $db_adaptor->get_GeneAdaptor(); # Bio::EnsEMBL::DBSQL::GeneAdaptor
my @genes = $gene_adaptor->fetch_all(); # Bio::EnsEMBL::Gene

#~ my $panH_genomeDB_adaptor = $registry->get_adaptor('pan_homology', 'compara', 'GenomeDB'); #Bio::EnsEMBL::Compara::DBSQL::GenomeDBAdaptor
#~ my $genome_db = $panH_genomeDB_adaptor->fetch_by_name_assembly("mus_musculus"); # Bio::EnsEMBL::Compara::GenomeDB
#~ my $db_adaptor = $genome_db->db_adaptor(); # Bio::EnsEMBL::DBSQL::DBAdaptor
#~ my $gene_adaptor = $db_adaptor->get_GeneAdaptor(); # Bio::EnsEMBL::DBSQL::GeneAdaptor
#~ my @genes = $gene_adaptor->fetch_all(); # Bio::EnsEMBL::Gene

#~ foreach my $gene (@genes){
	#~ foreach my $g (@{$gene}){
		#~ print $g->stable_id(),"\n";
	#~ }
#~ }

#~ foreach my $g (@{$genome_db}){
	#~ my $db_adaptor = $g->db_adaptor();
	#~ print $db_adaptor;
#~ }

#~ Bio::EnsEMBL::Registry->load_registry_from_db(
    #~ -host => 'mysql-eg-publicsql.ebi.ac.uk',
    #~ -user => 'anonymous',
    #~ -port => 5306);


#~ my $genome_db_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'GenomeDB');
#~ my $genome_db = $genome_db_adaptor->fetch_by_name_assembly("chlamydomonas_reinhardtii"); # Bio::EnsEMBL::Compara::GenomeDB
#~ my $db_adaptor = $genome_db->db_adaptor(); # Bio::EnsEMBL::DBSQL::DBAdaptor
#~ print $db_adaptor;
