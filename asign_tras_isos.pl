#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;
use common_stuff; #for build_lu
use gff3_annotation_stuff; #for parse_gff
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t 
This script is for idnetifying transcripts that came from
isoseq data. It uses exact coordinate matching to figure 
out which isoseq alignmets match the annotated transcripts

       Usage: asign_tras_isos.pl <annotations.gff> <alignments.gff>\n\n";

my $FILE_1 = $ARGV[0];
my $FILE_2 = $ARGV[1];
die($usage) unless $ARGV[0];
my %GENES;
my %ISO;

parse($FILE_1, 'g');
parse($FILE_2, 'i');
make_comparable_strings(\%GENES);
make_comparable_strings(\%ISO);
compare_strings();
#PostData(\%ISO);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub compare_strings{
    my %matches;
    
    foreach my $contig (keys %GENES){
	foreach my $tid (keys %{$GENES{$contig}}){
	    foreach my $aid (keys %{$ISO{$contig}}){
		next unless $GENES{$contig}{$tid}{'strand'} eq $ISO{$contig}{$aid}{'strand'};
		my $tcs = $GENES{$contig}{$tid}{'cs'};
		my $acs = $ISO{$contig}{$aid}{'cs'};
		if ($tcs eq $acs){
		    my $target = $ISO{$contig}{$aid}{'target'};
		    my $score  = $ISO{$contig}{$aid}{'score'};
		    push(@{$matches{$tid}}, [$aid, $score, $target]);
		}
	    }
	}
    }
    #dig through the matches if there are doubles look at the scores
    my %used; 
    foreach my $tid (keys %matches){
	my $num_matches = @{$matches{$tid}};
	my $isomatch = 'none';
	my $isoscore = 0;
	my $target = 'none';
	if ($num_matches > 1){
	    
            foreach my $pairs (@{$matches{$tid}}){
		if ($isoscore < $pairs->[1]){
		    $isomatch = $pairs->[0];
		    $isoscore = $pairs->[1];
		    $target   = $pairs->[2];
		}
		else {next;}
	    }
	}
	else {
	    $isomatch = $matches{$tid}->[0]->[0];
	    $isoscore = $matches{$tid}->[0]->[1];
	    $target   = $matches{$tid}->[0]->[2];
	}
	print "$tid\t$target\t$isomatch\t$isoscore\n" unless defined($used{$isomatch}); #exp
	$used{$isomatch} = 1;#exp
    }
}
#-----------------------------------------------------------------------------
sub make_comparable_strings{
    my $hr = shift;
    foreach my $contig (keys %$hr){
	foreach my $tid (keys %{$hr->{$contig}}){
	    #because of the multiparenting of exons they may not be 
	    #in the same order once they get into the arrays even though
	    #they have all of the same coordinates
	    my @sorted_ca = sort @{$hr->{$contig}->{$tid}->{'ca'}};
	    my $cs = join("_", @sorted_ca);
	    $hr->{$contig}->{$tid}->{'cs'} = $cs;
	}
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my $flag = shift;

    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /\#\#FASTA/;
	next if $line =~ /\#/;
	my @cols = split(/\t/, $line);
	if ($flag eq 'g'){
	    if ($cols[2] eq 'exon'){
		if ($cols[8] =~ /Parent=(\S+?);/ ||
		    $cols[8] =~ /Parent=(\S+)/){
		    my $parents = $1;
		    #remenber the multiple parenting of exons
		    my @tids = split(/\,/, $parents);
		    foreach my $tid (@tids){
			$GENES{$cols[0]}{$tid}{'strand'} = $cols[6];
			push(@{$GENES{$cols[0]}{$tid}{'ca'}}, ($cols[3], $cols[4]));
		    }
		}
	    }
	}
	
	elsif ($flag eq 'i'){
            if ($cols[2] eq 'match_part'){
                if ($cols[8] =~ /Parent=(\S+?);Target=(\S+)\s/){
                    my $aid = $1;
		    my $target = $2;
		    $ISO{$cols[0]}{$aid}{'target'} = $target;
		    $ISO{$cols[0]}{$aid}{'score'} = $cols[5];
		    $ISO{$cols[0]}{$aid}{'strand'} = $cols[6];
                    push(@{$ISO{$cols[0]}{$aid}{'ca'}}, ($cols[3], $cols[4]));
                }

            }

	}
	
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

