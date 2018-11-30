use Bio::EnsEMBL::Registry; #  import the Registry module which we use to establish this connection

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -port => 5306);

my $gene_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'GeneMember');
my $gene_member = $gene_member_adaptor->fetch_by_stable_id('ENSG00000004059');
#~ print $gene_member_adaptor; # Bio::EnsEMBL::Compara::DBSQL::GeneMemberAdaptor=HASH(0x1cd1ba8)
#~ print $gene_member; # Bio::EnsEMBL::Compara::GeneMember=HASH(0x1cd2118)
