use strict;
use warnings;
use Bio::EnsEMBL::Registry;
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
my $genome_db = $panH_genomeDB_adaptor->fetch_all_by_taxon_id("39947"); # Bio::EnsEMBL::Compara::GenomeDB
my $db_adaptor = $genome_db->db_adaptor();
