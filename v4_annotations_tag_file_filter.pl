#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
This script takes in the tag file that I created for the the v4 annotations
and assigns a confidence level. High confidence, medium confidence, and low
confidence.
\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];

my $data = parse($FILE);
my $ortholog_lists = populate_ortholog_list($data);
filter($data, $ortholog_lists);


#PostData($data);
report($data);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report{
    my $data_hr = shift;
    my $count=0;
    my $aed_count=0;
    foreach my $gid (keys %{$data_hr}){
	if (defined($data_hr->{$gid}->{'low_confidence'}) &&
	    $data_hr->{$gid}->{'low_confidence'} eq 'yes'){
	    $count++;
	}
	if (defined($data_hr->{$gid}->{'AED'}) &&
	    $data_hr->{$gid}->{'AED'} > 0.45){
	    $aed_count++;
	}
    }
#    print "aed > 0.45: $aed_count\n";
    print "filtered out: $count\n";
}
#-----------------------------------------------------------------------------
sub filter{
    my $data_hr = shift;
    my $o_lists_hr = shift;
    my $count=0;
    foreach my $gid (keys %{$data}){
	#I can change this line to see what can be rescued by orthology
	if (#$data_hr->{$gid}->{'covered_by_a_complete_TE'} eq 'yes' ||
	    #$data_hr->{$gid}->{'covered_by_a_fragmented_TE'} eq 'yes' ||
	    #$data_hr->{$gid}->{'AED'} > 0.45 ||
	    $data_hr->{$gid}->{'lenght_of_protein_product'} <= 50
	    ){
	    $count++;
	    if ($data_hr->{$gid}->{'brachypodium_distachyon_ortholog'} ne 'none'){
#from here
		if ($data_hr->{$gid}->{'brachypodium_distachyon_hom_typ'} =~ /ortholog_one2one/){
		    $data_hr->{$gid}->{'low_confidence'} = 'no';
		    $data_hr->{$gid}->{'rescued'} = 'only_ortholog';
		}
		else{
		    my @orths = split(/\,/, $data->{$gid}->{'brachypodium_distachyon_ortholog'});
		    foreach my $o (@orths){
			if ($o_lists_hr->{'brachypodium_distachyon_ortholog'}->{$o} == 1){ #important in subsequent iterations
			    $data_hr->{$gid}->{'low_confidence'} = 'no';
			    $data_hr->{$gid}->{'rescued'} = 'last_ortholog';
			} 
			else {$o_lists_hr->{'brachypodium_distachyon_ortholog'}->{$o}--;
			}
		    }
		}
#to here		
	    }
	    if ($data->{$gid}->{'oryza_sativa_ortholog'} ne 'none'){
#from here
		if ($data_hr->{$gid}->{'oryza_sativa_hom_typ'} =~ /ortholog_one2one/){
		    $data_hr->{$gid}->{'low_confidence'} = 'no';
		    $data_hr->{$gid}->{'rescued'} = 'only_ortholog';
		}
		else{
		    my @orths = split(/\,/, $data->{$gid}->{'oryza_sativa_ortholog'});
		    foreach my $o (@orths){
			if ($o_lists_hr->{'oryza_sativa_ortholog'}->{$o} == 1){ #important in subsequent iterations
			    $data_hr->{$gid}->{'low_confidence'} = 'no';
			    $data_hr->{$gid}->{'rescued'} = 'last_ortholog';
			} 
			else {$o_lists_hr->{'oryza_sativa_ortholog'}->{$o}--;
			}
		    }
		}
#to here				
		
	    }
	    if ($data_hr->{$gid}->{'setaria_italica_ortholog'} ne 'none'){
#from here
		if ($data_hr->{$gid}->{'setaria_italica_hom_typ'} =~ /ortholog_one2one/){
		    $data_hr->{$gid}->{'low_confidence'} = 'no';
		    $data_hr->{$gid}->{'rescued'} = 'only_ortholog';
		}
		else{
		    my @orths = split(/\,/, $data->{$gid}->{'setaria_italica_ortholog'});
		    foreach my $o (@orths){
			if ($o_lists_hr->{'setaria_italica_ortholog'}->{$o} == 1){ #important in subsequent iterations
			    $data_hr->{$gid}->{'low_confidence'} = 'no';
			    $data_hr->{$gid}->{'rescued'} = 'last_ortholog';
			} 
			else {$o_lists_hr->{'setaria_italica_ortholog'}->{$o}--;
			}
		    }
		}
#to here		
	    }
	    if ($data_hr->{$gid}->{'sorghum_bicolor_ortholog'} ne 'none'){
#from here
		if ($data_hr->{$gid}->{'sorghum_bicolor_hom_typ'} =~ /ortholog_one2one/){
		    $data_hr->{$gid}->{'low_confidence'} = 'no';
		    $data_hr->{$gid}->{'rescued'} = 'only_ortholog';
		}
		else{
		    my @orths = split(/\,/, $data->{$gid}->{'sorghum_bicolor_ortholog'});
		    foreach my $o (@orths){
			if ($o_lists_hr->{'sorghum_bicolor_ortholog'}->{$o} == 1){ #important in subsequent iterations
			    $data_hr->{$gid}->{'low_confidence'} = 'no';
			    $data_hr->{$gid}->{'rescued'} = 'last_ortholog';
			} 
			else {$o_lists_hr->{'sorghum_bicolor_ortholog'}->{$o}--;
			}
		    }
		}
#to here		
	    }
	    if ($data_hr->{$gid}->{'zea_mays_ortholog'} ne 'none'){
#from here
		if ($data_hr->{$gid}->{'zea_mays_hom_typ'} =~ /ortholog_one2one/){
		    $data_hr->{$gid}->{'low_confidence'} = 'no';
		    $data_hr->{$gid}->{'rescued'} = 'only_ortholog';
		}
		else{
		    my @orths = split(/\,/, $data->{$gid}->{'zea_mays_ortholog'});
		    foreach my $o (@orths){
			if ($o_lists_hr->{'zea_mays_ortholog'}->{$o} == 1){ #important in subsequent iterations
			    $data_hr->{$gid}->{'low_confidence'} = 'no';
			    $data_hr->{$gid}->{'rescued'} = 'last_ortholog';
			} 
			else {$o_lists_hr->{'zea_mays_ortholog'}->{$o}--;
			}
		    }
		}
#to here		
	    }
	    else {$data_hr->{$gid}{'low_confidence'} = 'yes' unless defined($data_hr->{$gid}{'low_confidence'});
	    }
	}
	#start here to look at other stuff
	
	
    
    }
print "total: $count\n";    
}
#-----------------------------------------------------------------------------
sub populate_ortholog_list{
    my $data_hr = shift;
    my %o_lists;

    foreach my $gid (keys %{$data}){
	foreach my $head (keys %{$data->{$gid}}){
	    if ($head =~ /ortholog$/){
		my @orths = split(/\,/, $data->{$gid}->{$head});
		foreach my $o (@orths){
		    $o_lists{$head}{$o}++;
		}
	    }
	}
    }
    return (\%o_lists);
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my @headers;
    my $num_cols;
    my %data;

    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	if ($line =~ /^gid/){
	    @headers = split(/\t/, $line);
	    $num_cols = @headers;
	}
	else {
	    my @cols = split(/\t/, $line);
	    
	    for (my $i = 1; $i < $num_cols; $i++){
		$data{$cols[0]}{$headers[$i]} = $cols[$i];
	    }
	}
    }
    $fh->close();
    return (\%data);
}
#-----------------------------------------------------------------------------

