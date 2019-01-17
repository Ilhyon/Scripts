use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $filename = '/home/anais/Documents/Data/Homology/Gene_all_sp.txt';
my @gene_list = ();

open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";
 
while (my $row = <$fh>) {
  chomp $row;
  push @gene_list, $row;
}

print "Read of the file done.\n";

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-eg-publicsql.ebi.ac.uk',
    -port => 4157,
    #~ -verbose => 1
);

print "Succeed to connect to the DB.\n";


# first you have to get a GeneMember object. In case of homology is a gene, in 
# case of family it can be a gene or a protein
my $pan_homology_adaptor  = $registry->get_adaptor('pan_homology', 'compara', 'GeneMember'); # Bio::EnsEMBL::Compara::DBSQL::GeneMemberAdaptor

# then you get the homologies where the member is involved

my $homology_adaptor = $registry->get_adaptor('pan_homology', 'compara', 'Homology'); #Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor=HASH(0x2148b100)

# That will return a reference to an array with all homologies (orthologues in
# other species and paralogues in the same one)
# Then for each homology, you can get all the Members implicated

my %Dico_homology;

foreach my $gene (@gene_list) {
	my $gene_member = $pan_homology_adaptor->fetch_by_stable_id($gene);
	my $homologies = $homology_adaptor->fetch_all_by_Member($gene_member);
	foreach my $homology (@{$homologies}) {
	  # You will find different kind of description
	  # see ensembl-compara/docs/docs/schema_doc.html for more details
	  my $pair_homologues = $homology->gene_list();
	  foreach my $gene (@{$pair_homologues}){
			my $id = $gene->stable_id();
			if($id ne 'ENSG00000100416' ){ # && $homology->description eq "ortholog_one2one"
				print $homology->description," ", $homology->taxonomy_level;
				print " ", $id, "\n";
			}
		}
	}
}






















