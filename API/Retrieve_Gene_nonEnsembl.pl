use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::Species;
use File::Path qw(make_path);
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-eg-publicsql.ebi.ac.uk',
    -port => 4157,
);

print "Succeed to connect to the DB\n";

sub get_genomes_by_name {
   my ($adaptor,@list) = @_;
   my @tmp = ();
   foreach my $sp (@list){
		my $genome_db = $adaptor->fetch_by_name_assembly($sp); # Bio::EnsEMBL::Compara::GenomeDB
		if($genome_db){ 
			push @tmp, $genome_db; 
		} else {
			print $sp,"\n";
		}
	}
   return @tmp;
}

sub get_genomes_by_taxonID {
   my ($adaptor,@list) = @_;
   my @tmp = ();
   foreach my $sp (@list){
		my $genome_db = $adaptor->fetch_all_by_taxon_id($sp); # Bio::EnsEMBL::Compara::GenomeDB
		if($genome_db){
			my $cpt = 0;
			foreach my $gene (@{$genome_db}){
				foreach my $g (@{$genome_db}){
					$cpt = $cpt +1;
					if(scalar @{$genome_db} == 1){
						push @tmp, $g;
					} else{
						print $g->name(),"\n";
					}
				}
			}
		} else {
			print $sp,"\n";
		}
	}
   return @tmp;
}

sub get_short_species_names {
	my ($sp) = @_;
	my $separator = "_";
	my $positionSecondWord = index($sp,$separator) + 1;
	my $initiales = substr($sp,0,1);
	$initiales = $initiales . substr($sp,$positionSecondWord,1);
	$initiales = uc($initiales);
	return $initiales;
}

sub check_file{
	my ($sp) = @_;
	my $directory = "~/Documents/Data/" . $sp;
	if(-e $directory and -d $directory){
		$directory = $directory . "/raw_material";
		unless(-e $directory and -d $directory){
			make_path($directory);
		}
	} else {
		make_path($directory);
		$directory = $directory . "/raw_material";
		make_path($directory);
		print $directory, "\n";
	}
}

sub write_output {
	my ($filename, @output_list) = @_;
	open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
	foreach my $line (@output_list){
		my $toPrint = $line ."\n";
		print $fh $toPrint;
	}
	close $fh;
}

sub get_gene_unspliced {
	my ($genome) = @_;
	my @tmp = ();
	my $db_adaptor = $genome->db_adaptor(); # Bio::EnsEMBL::DBSQL::DBAdaptor
	my $gene_adaptor = $db_adaptor->get_GeneAdaptor(); # Bio::EnsEMBL::DBSQL::GeneAdaptor
	my @genes = $gene_adaptor->fetch_all(); # Bio::EnsEMBL::Gene
	my $sp = $genome->name();
	my $initiales = get_short_species_names($sp);
	foreach my $gene (@genes){
		my $geneID = $gene->stable_id();
		my $strand = $gene->strand();
		my $start = $gene->start();
		my $end = $gene->end();
		my $header = $geneID . "|" . $strand . "|" . $start . "|" . $end . "\n";
		my $gene_unspliced = $gene->seq();
		my $fasta_gene_unspliced = $header . $gene_unspliced;
		push @tmp, $fasta_gene_unspliced; 
	}
	check_file($sp);
	my $filename = "~/Documents/Data/" . $sp . "/raw_material/" . $initiales . "_All_gene_unspliced.fasta";
	write_output($filename, @tmp);
}

sub get_all_genes {
	my ($genome) = @_;
	my $sp = $genome->name();
	my $initiales = get_short_species_names($sp);
	my $db_adaptor = $genome->db_adaptor(); # Bio::EnsEMBL::DBSQL::DBAdaptor
	my $gene_adaptor = $db_adaptor->get_GeneAdaptor(); # Bio::EnsEMBL::DBSQL::GeneAdaptor
	my @genes = $gene_adaptor->fetch_all(); # Bio::EnsEMBL::Gene
	check_file($sp);
	my $filename = "~/Documents/Data/" . $sp . "/raw_material/" . $initiales . "_All_GeneID.txt";
	write_output($filename, @genes);
}



my @list_sp_ensembl=("homo_sapiens","pan_troglodytes","pongo_abelii","mus_musculus","monodelphis_domestica","anolis_carolinensis","Ornithorhynchus_anatinus","gallus_gallus","danio_rerio","gasterosteus_aculeatus","xenopus_tropicalis","ciona_savignyi");
my @list_sp_metazoa=("caenorhabditis_elegans","drosophila_melanogaster","anopheles_gambiae","apis_mellifera","octopus_bimaculoides","pediculus_humanus","stegodyphus_mimosarum","amphimedon_queenslandica","brugia_malayi","mnemiopsis_leidyi","trichoplax_adhaerens");
my @list_sp_plants=("chlamydomonas_reinhardtii","chondrus_crispus","amborella_trichopoda","physcomitrella_patens","selaginella_moellendorffii","solanum_lycopersicum","vitis_vinifera","arabidopsis_thaliana");
my @list_sp_fungi=("aspergillus_nidulans","neurospora_crassa","saccharomyces_cerevisiae","schizosaccharomyces_pombe");
my @list_sp_protist=("bigelowiella_natans","dictyostelium_discoideum","emiliania_huxleyi","leishmania_major","tetrahymena_thermophila");
my @list_sp_bacteria=("francisella_tularensis_subsp_tularensis_schu_s4","escherichia_coli_str_k_12_substr_mg1655","bacillus_subtilis_subsp_subtilis_str_168","mycobacterium_tuberculosis_h37rv","enterococcus_faecalis_v583","mycoplasma_pneumoniae_m129","streptococcus_pneumoniae_TIGR4","borrelia_burgdorferi_B31","thermus_thermophilus_HB8","geobacter_sulfurreducens_PCA","wolbachia_endosymbiont_of_drosophila_melanogaster","aquifex_aeolicus_VF5","caulobacter_crescentus_CB15","cenarchaeum_symbiosum_a","coxiella_burnetii_rsa_493","gardnerella_vaginalis_0288E","helicobacter_pylori_26695","lactobacillus_plantarum_WCFS1","moraxella_catarrhalis_7169","myxococcus_xanthus_DK_1622","neisseria_meningitidis_Z2491","sulfolobus_solfataricus_P2","thermoplasma_acidophilum_DSM_1728","methanosarcina_acetivorans_C2A","pyrococcus_horikoshii_OT3","methanopyrus_kandleri_AV19","methanococcus_maripaludis_S2","aeropyrum_pernix_K1","archaeoglobus_fulgidus_DSM_4304","candidatus_korarchaeum_cryptofilum_OPF8","haloarcula_marismortui_ATCC_43049","halobacterium_salinarum_R1","hyperthermus_butylicus_DSM_5456");
my @taxonID_plants=("39947");
my @taxonID_bacteria=("158879","169963","272561","192222","262698","1131758","229193","212042","324602","716541","71421","272620","297246","203120","243277","228908","178306","420247");

my $genome_db_adaptor_metazoa = $registry->get_adaptor('metazoa', 'compara', 'GenomeDB');
my $genome_db_adaptor_plants = $registry->get_adaptor('plants', 'compara', 'GenomeDB');
my $genome_db_adaptor_fungi = $registry->get_adaptor('fungi', 'compara', 'GenomeDB');
my $genome_db_adaptor_protist = $registry->get_adaptor('protists', 'compara', 'GenomeDB');
my $genome_db_adaptor_bacteria = $registry->get_adaptor('bacteria', 'compara', 'GenomeDB');

my @genomes_db = ();

push @genomes_db, get_genomes_by_name($genome_db_adaptor_metazoa,@list_sp_metazoa); 
push @genomes_db, get_genomes_by_name($genome_db_adaptor_plants,@list_sp_plants); 
push @genomes_db, get_genomes_by_name($genome_db_adaptor_fungi,@list_sp_fungi); 
push @genomes_db, get_genomes_by_name($genome_db_adaptor_protist,@list_sp_protist); 
push @genomes_db, get_genomes_by_name($genome_db_adaptor_bacteria,@list_sp_bacteria);

push @genomes_db, get_genomes_by_taxonID($genome_db_adaptor_plants,@taxonID_plants); 
push @genomes_db, get_genomes_by_taxonID($genome_db_adaptor_bacteria,@taxonID_bacteria); 

#~ print scalar @genomes_db,"\n"; # 80

my @genomes_test = $genomes_db[1];

foreach my $genome (@genomes_test){
	get_all_genes($genome);
	get_gene_unspliced($genome);
}



































