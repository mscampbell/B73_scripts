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

\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];

parse($FILE);

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	next if $line =~ /^#/;
	my @cols = split(/\t/, $line);
	if ($cols[2] eq 'mRNA'){
	    $cols[1] = 'GMAP';
	    $cols[2] = 'expressed_sequence_match';
	    $cols[8] =~ s/Parent=\S+?;//;
	    $cols[8] =~ s/Parent=\S+?;//;
	    print join("\t", @cols), "\n";
	}
	elsif ($cols[2] eq 'exon'){
	    $cols[1] = 'GMAP';
	    $cols[2] = 'match_part';
	    print join("\t", @cols), "\n";
	}
	else {next;}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

