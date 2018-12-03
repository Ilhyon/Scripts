use strict;
use warnings;
use Bio::EnsEMBL::Registry;


my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-eg-publicsql.ebi.ac.uk',
    -port => 4157,
    #~ -verbose => 1
);
print "Connection to the data base Ensemble Genomes done.\n";

#~ print "Adapator done. \n";

#~ my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

#~ foreach my $db_adaptor (@db_adaptors) {
    #~ my $db_connection = $db_adaptor->dbc();

    #~ printf(
        #~ "species/group\t%s/%s\ndatabase\t%s\nhost:port\t%s:%s\n\n",
        #~ $db_adaptor->species(),   $db_adaptor->group(),
        #~ $db_connection->dbname(), $db_connection->host(),
        #~ $db_connection->port()
    #~ );
#~ }


# first you have to get a GeneMember object. In case of homology is a gene, in 
# case of family it can be a gene or a protein
my $pan_homology_adaptor  = $registry->get_adaptor('pan_homology', 'compara', 'GeneMember');
#~ print $pan_homology_adaptor # Bio::EnsEMBL::Compara::DBSQL::GeneMemberAdaptor
my $gene_member = $pan_homology_adaptor->fetch_by_stable_id('ENSG00000100416');
#~ print $gene_member; #Bio::EnsEMBL::Compara::GeneMember=HASH(0x1fff4598)

# then you get the homologies where the member is involved

my $homology_adaptor = $registry->get_adaptor('pan_homology', 'compara', 'Homology');
#~ print $homology_adaptor; #Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor=HASH(0x2148b100)
my $homologies = $homology_adaptor->fetch_all_by_Member($gene_member);
#~ print $homologies; #ARRAY(0x2144c178)

# That will return a reference to an array with all homologies (orthologues in
# other species and paralogues in the same one)
# Then for each homology, you can get all the Members implicated

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































#--------------- TEST -------------#

#~ my $specie = $pan_homology_adaptor->fetch_by_name_assembly('homo_sapiens');
#~ print $specie; # Bio::EnsEMBL::Compara::GenomeDB=HASH(0x21837440)

