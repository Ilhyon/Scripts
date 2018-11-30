use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-eg-publicsql.ebi.ac.uk',
    -port => 4157,
);

print "Connection to the data base Ensemble Genomes done.\n";

my $pan_homology_adaptor  = $registry->get_adaptor('pan_homology', 'compara', 'GeneMember');
my @gene_members = $pan_homology_adaptor->fetch_all();
#~ print @gene_members; ARRAY(0x1ac7b80)
#~ my $test = $pan_homology_adaptor->fetch_all();
#~ print $test; ARRAY(0xbeb7bf88)

foreach my $gene (@gene_members) {
  #~ print $gene; # array
  foreach my $i (@{$gene}){
	  print $i, "\n";
	  }
}

