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

my $cpt = 0;
print "0%";

$filename = '/home/anais/Documents/Data/Homology/Output.txt';
open($fh, '>', $filename) or die "Could not open file '$filename' $!";
foreach my $gene (@gene_list) {
	$cpt +=1;
	if($cpt == 64969){
		print "-10%";
	}elsif($cpt == 129939){
		print "--20%";
	}elsif($cpt == 194908){
		print "---30%";
	}elsif($cpt == 259877){
		print "----40%";
	}elsif($cpt == 324846){
		print "-----50%";
	}elsif($cpt == 389815){
		print "------60%";
	}elsif($cpt == 454784){
		print "-------70%";
	}elsif($cpt == 519753){
		print "--------80%";
	}elsif($cpt == 584722){
		print "---------90%";
	}elsif($cpt == 649696){
		print "----------100%";
	}
	my $gene_member = $pan_homology_adaptor->fetch_by_stable_id($gene);
	if($gene_member){
		my $homologies = $homology_adaptor->fetch_all_by_Member($gene_member);
		if($homologies){
			foreach my $homology (@{$homologies}) {
			  # You will find different kind of description
			  # see ensembl-compara/docs/docs/schema_doc.html for more details
			  my $pair_homologues = $homology->gene_list();
			  foreach my $genes (@{$pair_homologues}){
					my $id = $genes->stable_id();
					if($id ne $gene ){ # && $homology->description eq "ortholog_one2one"
						print $fh $gene, "\t", $homology->description,"\t", $id, "\n";
					}
				}
			}
		}
	}
}
 close $fh;





















