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
	my $directory = "/home/anais/Documents/Data/" . $sp;
	unless(-e $directory and -d $directory){
		make_path($directory);
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

sub get_transcript_info{
	my ($transcripts,$geneID,$strand) = @_;
	my @lines = ();
	foreach my $tr (@{$transcripts}){
		unless ((index($tr->seq_region_name(), "scaffold") != -1) && (index($tr->seq_region_name(), "contig") != -1)) {
			my $line = undef;
			my $five_UTR_start = undef;
			my $five_UTR_end = undef;
			my $three_UTR_start = undef;
			my $three_UTR_end = undef;
			my $trID = $tr->stable_id();
			my $CDS = $tr->get_all_CDS();
			my $exons = $tr->get_all_Exons();
			my $five_UTRs = $tr->get_all_five_prime_UTRs();
			my $three_UTRs = $tr->get_all_three_prime_UTRs();
			my $biotype = $tr->biotype();
			my $tr_start = $tr->start();
			my $tr_end = $tr->end();
			my $chr = $tr->seq_region_name();
			#~ print $chr,"\n";
			foreach my $five_UTR (@{$five_UTRs}){
				$five_UTR_start = $five_UTR->start();
				$five_UTR_end = $five_UTR->end();
			}
			foreach my $three_UTR (@{$three_UTRs}){
				$three_UTR_start = $three_UTR->start();
				$three_UTR_end = $three_UTR->end();
			}
			foreach my $exon (@{$exons}){
				my $exon_start = $exon->start();
				my $exon_end = $exon->end();
				my $rank = $exon->rank($tr);
				unless($five_UTR_start && $three_UTR_start){
					$line = $geneID . "\t" . $trID . "\t" . $biotype . "\t" . "\t" . "\t" . "\t" . "\t" . $exon_start . "\t" . $exon_end . "\t" . $rank . "\t" . $strand;
				}elsif($five_UTR_start && $three_UTR_start){
					$line = $geneID . "\t" . $trID . "\t" . $biotype . "\t" . $five_UTR_start . "\t" . $five_UTR_end . "\t" . $three_UTR_start . "\t" . $three_UTR_end . "\t" . $exon_start . "\t" . $exon_end . "\t" . $rank . "\t" . $strand;
				}elsif($five_UTR_start){
					$line = $geneID . "\t" . $trID . "\t" . $biotype . "\t" . $five_UTR_start . "\t" . $five_UTR_end . "\t" . "\t" . "\t" . $exon_start . "\t" . $exon_end . "\t" . $rank . "\t" . $strand;
				}elsif($three_UTR_start){
					$line = $geneID . "\t" . $trID . "\t" . $biotype . "\t" . "\t" . "\t" . $three_UTR_start . "\t" . $three_UTR_end . "\t" . $exon_start . "\t" . $exon_end . "\t" . $rank . "\t" . $strand;
				}
				push @lines, $line;
			}
		}
	}
	return @lines;
}

sub get_gene_unspliced {
	my ($genome) = @_;
	my @tmp = ();
	my @all_genes_ID = ();
	my @transcripts_info = ();
	my $sp = $genome->name();
	my $initiales = get_short_species_names($sp);
	my $db_adaptor = $genome->db_adaptor(); # Bio::EnsEMBL::DBSQL::DBAdaptor
	my $gene_adaptor = $db_adaptor->get_GeneAdaptor(); # Bio::EnsEMBL::DBSQL::GeneAdaptor
	my $genes = $gene_adaptor->fetch_all(); # Bio::EnsEMBL::Gene
	print $sp,"\n";
	my @transcript_unspliced = ();
	foreach my $gene (@{$genes}){
		my $geneID = $gene->stable_id();
		my $strand = $gene->strand();
		my $start = $gene->start();
		my $end = $gene->end();
		my $header = $geneID . "|" . $strand . "|" . $start . "|" . $end . "\n";
		my $gene_unspliced = $gene->seq();
		my $fasta_gene_unspliced = $header . $gene_unspliced;
		push @tmp, $fasta_gene_unspliced; 
		push @all_genes_ID, $geneID;
		my $transcripts = $gene->get_all_Transcripts();
		push @transcript_unspliced, get_transcript_info($transcripts,$geneID,$strand);
	}
	check_file($sp);
	my $filename = "/home/anais/Documents/Data/" . $sp . "/" . $initiales . "_gene_unspliced.txt";
	write_output($filename, @tmp);
	$filename = "/home/anais/Documents/Data/" . $sp . "/" . $initiales . "_GeneID.txt";
	write_output($filename, @all_genes_ID);
	$filename = "/home/anais/Documents/Data/" . $sp . "/" . $initiales . "_transcript_unspliced.txt";
	write_output($filename, @transcript_unspliced);
	
}

#~ my @list_sp_metazoa=("caenorhabditis_elegans","drosophila_melanogaster","anopheles_gambiae","apis_mellifera","octopus_bimaculoides","pediculus_humanus","stegodyphus_mimosarum","amphimedon_queenslandica","brugia_malayi","mnemiopsis_leidyi","trichoplax_adhaerens");
#~ my @list_sp_plants=("oryza_sativa","chlamydomonas_reinhardtii","chondrus_crispus","amborella_trichopoda","physcomitrella_patens","selaginella_moellendorffii","solanum_lycopersicum","vitis_vinifera","arabidopsis_thaliana");
#~ my @list_sp_fungi=("aspergillus_nidulans","neurospora_crassa","saccharomyces_cerevisiae","schizosaccharomyces_pombe");
#~ my @list_sp_protist=("bigelowiella_natans","dictyostelium_discoideum","emiliania_huxleyi","leishmania_major","tetrahymena_thermophila");
my @list_sp_bacteria=("methanococcus_maripaludis_S2","aeropyrum_pernix_K1","archaeoglobus_fulgidus_DSM_4304","candidatus_korarchaeum_cryptofilum_OPF8");
# "thermoplasma_acidophilum_DSM_1728","methanosarcina_acetivorans_C2A","campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819","brucella_abortus_bv_1_str_9_941","pseudomonas_aeruginosa_mpao1_p2","yersinia_pestis_biovar_microtus_str_91001","anaplasma_phagocytophilum_str_hz","chloroflexus_aurantiacus_j_10_fl","haemophilus_influenzae_rd_kw20","legionella_pneumophila_str_paris","vibrio_cholerae_o1_biovar_el_tor_str_n16961","pyrobaculum_aerophilum_str_im2","methanobrevibacter_smithii_atcc_35061","staphylococcus_aureus_subsp_aureus_n315","francisella_tularensis_subsp_tularensis_schu_s4","escherichia_coli_str_k_12_substr_mg1655","bacillus_subtilis_subsp_subtilis_str_168","mycobacterium_tuberculosis_h37rv","enterococcus_faecalis_v583","mycoplasma_pneumoniae_m129","borrelia_burgdorferi_B31","thermus_thermophilus_HB8","geobacter_sulfurreducens_PCA","wolbachia_endosymbiont_of_drosophila_melanogaster","aquifex_aeolicus_VF5","cenarchaeum_symbiosum_a","gardnerella_vaginalis_0288E","myxococcus_xanthus_DK_1622","neisseria_meningitidis_Z2491",
my @error = ("chlamydia_trachomatis_d_uw_3_cx","nanoarchaeum_equitans_kin4_m","streptococcus_pneumoniae_TIGR4","sulfolobus_solfataricus_P2","pyrococcus_horikoshii_OT3","halobacterium_salinarum_R1");

#~ my $genome_db_adaptor_metazoa = $registry->get_adaptor('metazoa', 'compara', 'GenomeDB');
#~ my $genome_db_adaptor_plants = $registry->get_adaptor('plants', 'compara', 'GenomeDB');
#~ my $genome_db_adaptor_fungi = $registry->get_adaptor('fungi', 'compara', 'GenomeDB');
#~ my $genome_db_adaptor_protist = $registry->get_adaptor('protists', 'compara', 'GenomeDB');
my $genome_db_adaptor_bacteria = $registry->get_adaptor('bacteria', 'compara', 'GenomeDB');

my @genomes_db = ();

#~ push @genomes_db, get_genomes_by_name($genome_db_adaptor_metazoa,@list_sp_metazoa); 
#~ push @genomes_db, get_genomes_by_name($genome_db_adaptor_plants,@list_sp_plants); 
#~ push @genomes_db, get_genomes_by_name($genome_db_adaptor_fungi,@list_sp_fungi); 
#~ push @genomes_db, get_genomes_by_name($genome_db_adaptor_protist,@list_sp_protist); 
push @genomes_db, get_genomes_by_name($genome_db_adaptor_bacteria,@list_sp_bacteria);

#~ print scalar @genomes_db,"\n"; # 80

#~ my @genomes_test = $genomes_db[1];
#print $genomes_db[1],"\n";

#~ my $test = $genome_db_adaptor_bacteria->fetch_by_name_assembly("campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819"); # Bio::EnsEMBL::Compara::GenomeDB
		
#~ get_gene_unspliced($test);
foreach my $genome (@genomes_db){
	get_gene_unspliced($genome);
}



































