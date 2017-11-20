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
my %DATA;
my %mRNA_LU;
parse($FILE);
report();
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report{
    foreach my $key (keys %DATA){
	if ($DATA{$key} == 1){
	    print $key,"\n" if $mRNA_LU{$key};
	    
	}
    }

}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /\#\#FASTA/;
	next if $line =~ /\#/;
	my @cols = split(/\t/, $line);
	if ($cols[2] eq 'exon'){
	    if  ($line =~ /Parent=(\S+)/){
		my @par = split(/,/, $1);
		foreach my $x (@par){
		    $DATA{$x}++;
		}
	    }
	}
	if ($cols[2] eq 'mRNA'){
	    my ($id) = $line =~ /ID=(\S+?);/;
	    $mRNA_LU{$id}++;
	}
	
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

