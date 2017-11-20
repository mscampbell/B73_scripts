#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_r $opt_t $opt_g $opt_p $opt_c $opt_f $opt_l $opt_o $opt_s $opt_b $opt_i $opt_z);
getopts('r:t:g:p:c:f:l:o:s:b:i:z:');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
This script is for interogating the v4 gene models. I am planning to parse
the QI portions of the maker output and looks at the completeness of the gene 
models. I'm also going to look at the support for the single exon genes. 

explore_v4_annotatons.pl -g test_gene.gff -r test_reference.fasta

\n\n";

my $GFF    = $opt_g;
my $REF    = $opt_r;
my $R_TRAN = $opt_t;
my $TE_PF  = $opt_p;
my $TE_C   = $opt_c;
my $TE_F   = $opt_f;
my $LNC    = $opt_l;
my $OS_O   = $opt_o; #rice orthology
my $SB_O   = $opt_s; #sorghum orthology
my $BD_O   = $opt_b;
my $SI_O   = $opt_i;
my $V3_O   = $opt_z;
die($usage) unless ($opt_g);
my %BIG_HASH;

my $gff_hr = parse_gff($GFF);
my $ref_fasta_hr = parse_fasta($REF) if $opt_r;
my ($incomplete_cds_ar) = identify_incomplete_transcripts($gff_hr, $ref_fasta_hr);
my $single_ex_genes_ar = identify_low_quality_single_exon_genes($gff_hr);
genes_with_pfam_blast($gff_hr, $TE_PF); #identify the TE domains in here
genes_covered_by_TE_or_lnc($TE_C, 'c');
genes_covered_by_TE_or_lnc($TE_F, 'f');
genes_covered_by_TE_or_lnc($LNC, 'l');
rep_transcript_data($gff_hr, $R_TRAN);
add_orthology($OS_O);
add_orthology($SB_O);
add_orthology($BD_O);
add_orthology($SI_O);
add_orthology($V3_O);

#report_incomplete_models($incomplete_cds_ar);
make_the_tag_file(\%BIG_HASH);

#PostData($gff_hr);
#check_stuff();
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub make_the_tag_file{
# these are potential keys in the big hash not all will be present on every
# gene
#############GENE CENTRIC INFO#############
#AED
#start_with_no_stop
#stop_with_no_start
#complete_cds
#no_start_or_stop
#frac_ss_confirmed_by_est
#frac_exons_with_est_or_pro_overlap
#############SPPROT AND PFAM#############
#similarity_to_sprot
#contains_non_TE_pfam_domain
#############TE#############
#contains_TE_pfam_domain
#covered_by_a_complete_TE
#covered_by_a_fragmented_TE
#############lncRNA#############
#covered_by_a_lncRNA
#############orthology#############
#brachypodium_distachyon_ortholog
#brachypodium_distachyon_hom_typ
#oryza_sativa_ortholog
#oryza_sativa_hom_typ
#setaria_italica_ortholog
#setaria_italica_hom_typ
#sorghum_bicolor_ortholog
#sorghum_bicolor_hom_typ
#zea_mays_ortholog
#zea_mays_hom_typ

    
    print "gid\tAED\tCDS_status\tfrac_ss_confirmed_by_est\tlenght_of_protein_product\tfrac_exons_with_est_or_pro_overlap\tsimilarity_to_sprot\tcontains_non_TE_pfam_domain\tcontains_TE_pfam_domain\tcovered_by_a_complete_TE\tcovered_by_a_fragmented_TE\tcovered_by_a_lncRNA\tbrachypodium_distachyon_ortholog\tbrachypodium_distachyon_hom_typ\toryza_sativa_ortholog\toryza_sativa_hom_typ\tsetaria_italica_ortholog\tsetaria_italica_hom_typ\tsorghum_bicolor_ortholog\tsorghum_bicolor_hom_typ\tzea_mays_ortholog\tzea_mays_hom_typ\n";
    my $hr = shift;
    foreach my $gid (keys %{$hr}){
	print "$gid\t"; #gid
	print "$hr->{$gid}->{'AED'}\t"; #AED
	#CDS_status
	if ($hr->{$gid}->{'complete_cds'}){
	    print "complete_cds\t";
	}
	elsif ($hr->{$gid}->{'stop_with_no_start'}){
	    print "stop_with_no_start\t"
	}
	elsif ($hr->{$gid}->{'start_with_no_stop'}){
	    print "start_with_no_stop\t"
	}
	elsif ($hr->{$gid}->{'no_start_or_stop'}){
	    print "no_start_or_stop\t"
	}
	else {die "you missed a cds posibility\n";}

	#frac_ss_confirmed_by_est
	print "$hr->{$gid}->{'frac_ss_confirmed_by_est'}\t"; 
	
	#lenght_of_protein_product
	print "$hr->{$gid}->{'lenght_of_protein_product'}\t"; 

        #frac_exons_with_est_or_pro_overlap
	print "$hr->{$gid}->{'frac_exons_with_est_or_pro_overlap'}\t";

	#similarity_to_sprot
	print "$hr->{$gid}->{'similarity_to_sprot'}\t"; 

	#contains_non_TE_pfam_domain
	if(defined($hr->{$gid}->{'contains_non_TE_pfam_domain'})){
	    print "$hr->{$gid}->{'contains_non_TE_pfam_domain'}\t"; 
	}
	else{print "no\t";}
	
	#contains_TE_pfam_domain
	if (defined($hr->{$gid}->{'contains_TE_pfam_domain'})){
	    print "$hr->{$gid}->{'contains_TE_pfam_domain'}\t"; 
	}
	else{print "no\t";}
	
        #covered_by_a_complete_TE
	if (defined($hr->{$gid}->{'covered_by_a_complete_TE'})){
	    print "$hr->{$gid}->{'covered_by_a_complete_TE'}\t"; 
	}
	else{print "no\t";}

        #covered_by_a_fragmented_TE
	if(defined($hr->{$gid}->{'covered_by_a_fragmented_TE'})){
	    print "$hr->{$gid}->{'covered_by_a_fragmented_TE'}\t"; 
	}
	else{print "no\t";}
	
        #covered_by_a_lncRNA
	if (defined($hr->{$gid}->{'covered_by_a_lncRNA'})){
	    print "$hr->{$gid}->{'covered_by_a_lncRNA'}\t"; 
	}
	else{print "no\t";}
       
        #brachypodium_distachyon_ortholog
        #brachypodium_distachyon_hom_typ
	if (defined($hr->{$gid}->{'brachypodium_distachyon_ortholog'})){
	    print join(",", @{$hr->{$gid}->{'brachypodium_distachyon_ortholog'}}), "\t";
	    print join(",", @{$hr->{$gid}->{'brachypodium_distachyon_hom_typ'}}), "\t";
	}
	else{print "none\tnone\t";}
	
        #oryza_sativa_ortholog
        #oryza_sativa_hom_typ
	if (defined($hr->{$gid}->{'oryza_sativa_ortholog'})){
	    print join(",", @{$hr->{$gid}->{'oryza_sativa_ortholog'}}), "\t";
	    print join(",", @{$hr->{$gid}->{'oryza_sativa_hom_typ'}}), "\t";
	}
	else{print "none\tnone\t";}

        #setaria_italica_ortholog
        #setaria_italica_hom_typ
	if (defined($hr->{$gid}->{'setaria_italica_ortholog'})){
	    print join(",", @{$hr->{$gid}->{'setaria_italica_ortholog'}}), "\t";
	    print join(",", @{$hr->{$gid}->{'setaria_italica_hom_typ'}}), "\t";
	}
	else{print "none\tnone\t";}
        
        #sorghum_bicolor_ortholog
        #sorghum_bicolor_hom_typ
	if (defined($hr->{$gid}->{'sorghum_bicolor_ortholog'})){
	    print join(",", @{$hr->{$gid}->{'sorghum_bicolor_ortholog'}}), "\t";
	    print join(",", @{$hr->{$gid}->{'sorghum_bicolor_hom_typ'}}), "\t";
	}
	else{print "none\tnone\t";}

        #zea_mays_ortholog
        #zea_mays_hom_typ
	if (defined($hr->{$gid}->{'zea_mays_ortholog'})){
	    print join(",", @{$hr->{$gid}->{'zea_mays_ortholog'}}), "\t";
	    print join(",", @{$hr->{$gid}->{'zea_mays_hom_typ'}});
	}
	else{print "none\tnone";}

	print "\n";
    }
}
#-----------------------------------------------------------------------------
sub add_orthology{
    my $file = shift;
    my $fh = new FileHandle;

    $fh->open($file);
    while (defined(my $line = <$fh>)){
        chomp($line);
	next if $line =~ /^gid/;
	my @cols = split(/\t/, $line);
	my $gid = $cols[1];
	my $oid = $cols[0];
	my $o_relationship = $cols[8];
	my $species = $cols[4];
	#$species .= '_ortholog';
	
	push(@{$BIG_HASH{$gid}{$species.'_ortholog'}},  $oid);
	push(@{$BIG_HASH{$gid}{$species.'_hom_typ'}}, $o_relationship);
    }
}
#-----------------------------------------------------------------------------
sub rep_transcript_data{
    my $gff_hr      = shift;
    my $rep_t_file  = shift;
    my %rep_t_lu;

    my $fh = new FileHandle;
    $fh->open($rep_t_file);
    while (defined(my $line = <$fh>)){
        chomp($line);
	$rep_t_lu{$line}=1;
    }
    
    foreach my $gid(keys %{$gff_hr}){
        foreach my $tid  (keys %{$gff_hr->{$gid}->{'tl'}}){
	    if (defined($rep_t_lu{$tid})){
		my $aed = $gff_hr->{$gid}->{'tl'}->{$tid}->{'_AED'};
		my $f_ss_c = $gff_hr->{$gid}->{'tl'}->{$tid}->{'frac_ss_confirmed_by_est'};
		my $fp_utr = $gff_hr->{$gid}->{'tl'}->{$tid}->{'five_prime_utr_length'};
		my $tp_utr = $gff_hr->{$gid}->{'tl'}->{$tid}->{'three_prime_utr_length'};
		my $ex_with_sup = $gff_hr->{$gid}->{'tl'}->{$tid}->{'frac_exons_with_est_or_pro_overlap'};	
		my $tran_cds_lenght = $gff_hr->{$gid}->{'tl'}->{$tid}->{'lenght_of_protein_product'};

		$BIG_HASH{$gid}{'AED'} = $aed;
		$BIG_HASH{$gid}{'frac_ss_confirmed_by_est'} = $f_ss_c;
		$BIG_HASH{$gid}{'frac_exons_with_est_or_pro_overlap'} = $ex_with_sup;
		$BIG_HASH{$gid}{'lenght_of_protein_product'} = $tran_cds_lenght;
	    }
	}
    }
}
#-----------------------------------------------------------------------------
sub genes_covered_by_TE_or_lnc{
    my $id_file = shift;
    my $flag = shift;

    my $fh = new FileHandle;
    $fh->open($id_file);
    while (defined(my $line = <$fh>)){
        chomp($line);
	my @cols = split(/\t/, $line);
	if ($flag eq 'c'){
	    $BIG_HASH{$cols[0]}{'covered_by_a_complete_TE'}= 'yes';
	}
	elsif ($flag eq 'f'){
	    unless (defined($BIG_HASH{$cols[0]}{'covered_by_a_complete_TE'}) &&
		    $BIG_HASH{$cols[0]}{'covered_by_a_complete_TE'} eq 'yes'){
		$BIG_HASH{$cols[0]}{'covered_by_a_fragmented_TE'}= 'yes';
	    }
	}
	elsif ($flag = 'l'){
	    	$BIG_HASH{$cols[0]}{'covered_by_a_lncRNA'}= 'yes';
	}
    }
    
    $fh->close();
    
	
}
#-----------------------------------------------------------------------------
sub genes_with_pfam_blast{
    my $gff_hr    = shift;
    my $te_d_file = shift;
    my %te_d_lu;
    my $fh = new FileHandle;
    $fh->open($te_d_file);
    
    while (defined(my $line = <$fh>)){
        chomp($line);
	$te_d_lu{$line}=1;
    }
    $fh->close();
    
    foreach my $gid(keys %{$gff_hr}){
	#get the blast data
	if ($gff_hr->{$gid}->{'gl'}->{'Note'} =~ /^Similar/){
	    $BIG_HASH{$gid}{'similarity_to_sprot'} = 'yes';
	}
	else {
	    $BIG_HASH{$gid}{'similarity_to_sprot'} = 'no';
	}
	#get the TE domain data
	foreach my $tid  (keys %{$gff_hr->{$gid}->{'tl'}}){
	    if (defined($gff_hr->{$gid}->{'tl'}->{$tid}->{'Dbxref'})){
		my @pfams = split(/\,/, $gff_hr->{$gid}->{'tl'}->{$tid}->{'Dbxref'});
		foreach my $pfid (@pfams){
		    $pfid =~ s/Pfam://;
		    if (defined($te_d_lu{$pfid})){
			$BIG_HASH{$gid}{'contains_TE_pfam_domain'}++;
		}
		    else {
			$BIG_HASH{$gid}{'contains_non_TE_pfam_domain'}++;
		    }
		}
	    }
	    #could add another conditonal here and go into the QI stats or AED
	}
    }    
}
#-----------------------------------------------------------------------------
sub check_stuff{
    my @start_stop;
    my $strand        = '-';
    my $ref_seq       = 'AAAACCCCAAAATTTTGGGG';
                        #12345678901234567890 
    push(@start_stop, [3,9]); 
    push(@start_stop, [12,17]);
    push(@start_stop, [19,20]);

    #AACCCCAATTTTGGG
    my @sorted_start_stop = sort {$a->[0] <=> $b->[0]} @start_stop;
    my $cds = get_discontinuous_feature_sub_seq(\@sorted_start_stop, $strand, $ref_seq);
    #print "$cds\n";
}
#-----------------------------------------------------------------------------
sub report_incomplete_models{
    my $incomplete_cds_ar = shift;
    
    foreach my $ar (@{$incomplete_cds_ar}){
	print STDOUT join("\t", @{$ar}), "\n" if (defined(@{$ar}));
    }
    
}
#-----------------------------------------------------------------------------
sub get_discontinuous_feature_sub_seq{
    my $start_stop_ar = shift;
    my $strand        = shift;
    my $ref_seq       = shift;
    my $subfeature_seq_pos;
    my $subfeature_seq;

    foreach my $beg_end (@{$start_stop_ar}){
	my $begin = $beg_end->[0];
	my $end = $beg_end->[1];
	my $length =  $end - $begin +1;
	my $index  = $begin -1;
	$subfeature_seq_pos .= substr($ref_seq, $index, $length);
#	print "$begin\t$end\n";
    }
    if ($strand eq '-'){
	$subfeature_seq = reverse_comp($subfeature_seq_pos);
    }
    elsif ($strand eq '+'){
	$subfeature_seq = $subfeature_seq_pos;
    }
 #   print "$subfeature_seq_pos\n$subfeature_seq\n\n";
    return ($subfeature_seq);
}
#-----------------------------------------------------------------------------
sub identify_incomplete_transcripts{
    my $gff_hr = shift;
    my $ref_fasta_hr = shift;
    my @incomplet_cds_tids;
    foreach my $gid (keys %$gff_hr){
	    #get the reference sequence id
	    my $ref_seq_id = $gff_hr->{$gid}->{'gl'}->{'refseq'};
	    #use the id to get the actual sequence
	    my $ref_seq = $ref_fasta_hr->{$ref_seq_id};
	    my $strand = $gff_hr->{$gid}->{'gl'}->{'strand'};
	foreach my $tid (keys %{$gff_hr->{$gid}{'tl'}}){
	    my @start_stop;
	    my $begin_end_phase_ar = $gff_hr->{$gid}->{'tl'}->{$tid}->{'CDS'}->{'coordinates_and_pahse'};

	    foreach my $array_ref (@{$begin_end_phase_ar}){
		push(@start_stop, [$array_ref->[0], $array_ref->[1]]); 
	    }
	    #print "$gid is bad\n" unless defined($ref_seq);
#	    print "$tid\n";
	    #sort the start stop array (this may have bit me)
	    my @sorted_start_stop = sort {$a->[0] <=> $b->[0]} @start_stop;
	    my $cds = get_discontinuous_feature_sub_seq(\@sorted_start_stop, $strand, $ref_seq);
	    if ($cds =~ /^ATG/ && 
		($cds =~ /TAG$/ || $cds =~ /TAA$/ || $cds =~ /TGA$/)){
		#has a connoical start and stop
		push(@incomplet_cds_tids, [$tid, $gid, $strand, $ref_seq_id, 'complete']);
#		next;
		$BIG_HASH{$gid}{'complete_cds'}++;
	    }
	    elsif ($cds =~/^ATG/ && $cds !~ /TAG$/ && $cds !~ /TAA$/ && $cds !~ /TGA$/){
	    #has a start but no stop
		push(@incomplet_cds_tids, [$tid, $gid, $strand, $ref_seq_id, 'start_with_no_stop']);
		$BIG_HASH{$gid}{'start_with_no_stop'}++;
	    }
	    elsif ($cds !~/^ATG/ &&
		   ($cds =~ /TAG$/ || $cds =~ /TAA$/ || $cds =~ /TGA$/)){
		#Doesn't have a start but has a stop
		push(@incomplet_cds_tids, [$tid, $gid, $strand, $ref_seq_id, 'stop_with_no_start']);
		$BIG_HASH{$gid}{'stop_with_no_start'}++;
	    }
	    elsif ($cds !~/^ATG/ && $cds !~ /TAG$/ && $cds !~ /TAA$/ && $cds !~ /TGA$/){
		#Doesn't have a start or a stop
		push(@incomplet_cds_tids, [$tid, $gid, $strand, $ref_seq_id, 'no_start_or_stop']);
		$BIG_HASH{$gid}{'no_start_or_stop'}++;
	    }
	    else {die "you missed a possible start and stop combination";}
	}
    }
    return (\@incomplet_cds_tids);
}
#-----------------------------------------------------------------------------
sub identify_low_quality_single_exon_genes{
    my $gff_hr = shift;
    my @single_ex_genes;
    foreach my $gid (keys %$gff_hr){
	my $total_trans =0;
	my $single_exon_trans =0;
	foreach my $tid (keys %{$gff_hr->{$gid}{'tl'}}){
	    $total_trans++;
	    my $num_exons = $gff_hr->{$gid}->{'tl'}->{$tid}->{'number_of_exons'};
	    if ($num_exons == 1){
		$single_exon_trans++;
	    } 
	}
	#print "$total_trans == $single_exon_trans\n";
	if ($total_trans == $single_exon_trans){
	    push(@single_ex_genes, $gid);
	}
    }
    return (\@single_ex_genes);
}
#-----------------------------------------------------------------------------
sub parse_fasta{
    print STDERR "Now parsing Fasta sequence\n";
    my $file = shift;
    my %ref_fasta;
    my $def;
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	
	if ($line =~ /^>(\S+)\s/ || $line =~ /^>(\S+)/){
	    $def = $1;
	} 

	else{
	    $ref_fasta{$def} .= $line;
	}
    }
    $fh->close;
    return (\%ref_fasta);
}
#-----------------------------------------------------------------------------
sub parse_gff{
    print STDERR "Now Parsing GFF\n";
    my $file = shift;
    my %gff;
    my %gid_tid_lu;
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	my @cols = split(/\t/ , $line);

#parse the gene lines
	if ($cols[2] eq 'gene'){
	    my $gid;
	    if ($cols[8] =~ /ID=(\S+?);/ || $cols[8] =~ /ID=(\S+)/){
		$gid = $1;
	    }
	    else {die "Not valid gff3";}
	    #load everything before column 9
	    $gff{$gid}{'gl'}{'begin'}   = $cols[3];
	    $gff{$gid}{'gl'}{'end'}     = $cols[4];
	    $gff{$gid}{'gl'}{'strand'}  = $cols[6];
	    $gff{$gid}{'gl'}{'refseq'}  = $cols[0];
	    $gff{$gid}{'gl'}{'source'}  = $cols[1];
	    $gff{$gid}{'gl'}{'score'}   = $cols[5];
	    $gff{$gid}{'gl'}{'phase'}   = $cols[7];
	    
            #load column 9
	    my @col9 = split(/;/, $cols[8]);
	    foreach my $kvp (@col9){
		my ($k, $v) = split(/=/, $kvp);
		$gff{$gid}{'gl'}{$k} = $v;
	    }
	}
#parse the mRNA lines
	elsif ($cols[2] eq 'mRNA'){
	    my $tid;
	    my $gid;
	    if ($cols[8] =~ /ID=(\S+?);/){
		$tid = $1;
	    }
	    else {die "Not valid gff3";}
	    
	    if ($cols[8] =~ /Parent=(\S+?);/ || $cols[8] =~ /Parent=(\S+)/){
		$gid = $1;
	    }
	    else {die "Not valid gff3";}
	    #build a look up for later
	    $gid_tid_lu{'g2t'}{$gid} = $tid;
	    $gid_tid_lu{'t2g'}{$tid} = $gid;

            #load everything before column 9
	    $gff{$gid}{'tl'}{$tid}{'begin'}   = $cols[3];
	    $gff{$gid}{'tl'}{$tid}{'end'}     = $cols[4];
	    $gff{$gid}{'tl'}{$tid}{'strand'}  = $cols[6];
	    $gff{$gid}{'tl'}{$tid}{'refseq'}  = $cols[0];
	    $gff{$gid}{'tl'}{$tid}{'source'}  = $cols[1];
	    $gff{$gid}{'tl'}{$tid}{'score'}   = $cols[5];
	    #$gff{$gid}{'tl'}{$tid}{'phase'}   = $cols[7];
	    
            #load column 9
	    my @col9 = split(/;/, $cols[8]);
	    foreach my $kvp (@col9){
		my ($k, $v) = split(/=/, $kvp);
		
		if ($k eq '_QI'){
		    my @qi = split(/\|/, $v);

		    $gff{$gid}{'tl'}{$tid}{'five_prime_utr_length'} = $qi[0]; 
		    $gff{$gid}{'tl'}{$tid}{'frac_ss_confirmed_by_est'} = $qi[1]; 
		    $gff{$gid}{'tl'}{$tid}{'frac_exons_matching_est'} = $qi[2]; 
		    $gff{$gid}{'tl'}{$tid}{'frac_exons_with_est_or_pro_overlap'} = $qi[3]; 
		    $gff{$gid}{'tl'}{$tid}{'frac_ss_confirmed_by_gene_pred'} = $qi[4]; 
		    $gff{$gid}{'tl'}{$tid}{'frac_exons_with_gene_pred_overlap'} = $qi[5]; 
		    $gff{$gid}{'tl'}{$tid}{'number_of_exons'} = $qi[6]; 
		    $gff{$gid}{'tl'}{$tid}{'three_prime_utr_length'} = $qi[7]; 
		    $gff{$gid}{'tl'}{$tid}{'lenght_of_protein_product'} = $qi[8]; 

		}
		else{
		    $gff{$gid}{'tl'}{$tid}{$k} = $v;
		}
	    }
	    
	}
#parse the exon, cds, and utr lines
	elsif ($cols[2] eq 'exon' ||
	       $cols[2] eq 'CDS' ||
	       $cols[2] eq 'three_prime_UTR' ||
	       $cols[2] eq 'five_prime_UTR'){
	       
	    my @tids; #becasue of multiple parenting this my have more
                     #than one id in it (comma separated). If not this
	             #will jsut turn into an array with one value
	    if ($cols[8] =~ /Parent=(\S+?);/ || $cols[8] =~ /Parent=(\S+)/){
                @tids = split(/,/, $1);

            }
            else {die "Not valid gff3";}
	    #get the featrue ids of things that are not exons
	    my $fid;
	    if (($cols[2] eq 'CDS' ||
		 $cols[2] eq 'three_prime_UTR' ||
		 $cols[2] eq 'five_prime_UTR') && $cols[8] =~ /ID=(\S+?);/){
                $fid = $1;
            }
	    
            else {die " I don't think this is valid gff3. I need ids for cds and utr features\n$cols[2] $cols[8]\n" unless $cols[2] eq 'exon';}
	    
	    foreach my $tid (@tids){
		my $gid = $gid_tid_lu{'t2g'}{$tid};
		if ($cols[2] eq 'exon'){
		    push (@{$gff{$gid}{'tl'}{$tid}{'exons'}}, [$cols[3], $cols[4]]);

		}
		else {
		    $gff{$gid}{'tl'}{$tid}{$cols[2]}{'id'} = $fid;
		    if ($cols[2] eq 'CDS'){
			push (@{$gff{$gid}{'tl'}{$tid}{$cols[2]}{'coordinates_and_pahse'}}, [$cols[3], $cols[4], $cols[7]]);
		    }
		    else{
			push (@{$gff{$gid}{'tl'}{$tid}{$cols[2]}{'coordinates'}}, [$cols[3], $cols[4]]);
		    }
		}
	    }
	}
	
	else {
	    die "unexpected typ: $cols[2]\n";
	}
    }
    $fh->close();
    return (\%gff);
}
#-----------------------------------------------------------------------------
sub reverse_comp{

    my $seq = shift;
    my $r_seq;

    $r_seq = reverse($seq);
    $r_seq =~ tr/ACGTYRKMB/TGCARUMKV/;

    return $r_seq;
}
#-----------------------------------------------------------------------------
