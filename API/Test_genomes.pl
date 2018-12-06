use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::Species;
#~ use Bio::EnsEMBL::LookUp;
#~ my $lookup = Bio::EnsEMBL::LookUp->new();
#~ my $dba = $lookup->get_by_name_exact('escherichia_coli_str_k_12_substr_mg1655');   
#~ my @dbas = @{$lookup->get_all_by_name_pattern('escherichia_coli_.*')};
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-eg-publicsql.ebi.ac.uk',
    -port => 4157,
);

my $panH_genomeDB_adaptor = $registry->get_adaptor('pan_homology', 'compara', 'GenomeDB');
my $genome_db = $panH_genomeDB_adaptor->fetch_by_name_assembly("haloarcula_marismortui_ATCC_43049"); # Bio::EnsEMBL::Compara::GenomeDB
my $db_adaptor = $genome_db->db_adaptor(); # Bio::EnsEMBL::DBSQL::DBAdaptor
my $gene_adaptor = $db_adaptor->get_GeneAdaptor(); # Bio::EnsEMBL::DBSQL::GeneAdaptor
my @genes = $gene_adaptor->fetch_all(); # Bio::EnsEMBL::Gene

foreach my $gene (@genes){
	foreach my $g (@{$gene}){
		print $g->stable_id(),"\n";
	}
}
