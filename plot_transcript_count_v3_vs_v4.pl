#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_s $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('segpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t

\n\n";

my $GFF1 = $ARGV[0];
my $GFF2 = $ARGV[1];
my $MAP  = $ARGV[2];

die($usage) unless $ARGV[0];
my %DATA;
parse_gff($GFF1);
parse_gff($GFF2);
count_and_print($MAP);
#PostData(\%DATA);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub count_and_print{
    my $file = shift;
    my $v4_has_more = 0;
    my $no_change   = 0;
    my $v4_has_less = 0;
    my $total       = 0;
    print "v4_count\tv3_count\n";
    my $fh = new FileHandle;
    $fh->open($file);

    while (defined(my $line = <$fh>)){
        chomp($line);
	next if $line =~ /^\#/;
	$total++;

	my @cols = split(/\t/, $line);
	my $v3id = $cols[0];
	my $v4id = $cols[2];
	#print "$v3id\t$v4id\n";
	next unless (defined($DATA{$v3id}) && defined($DATA{$v4id}));
	my $v3count = $DATA{$v3id};
	my $v4count = $DATA{$v4id};
	
	$v4_has_more++ if $v4count > $v3count;
	$no_change++ if $v4count == $v3count;
	$v4_has_less++ if $v4count < $v3count;
	print "$v4count\t$v3count\n";
    }
    $fh->close();
    if ($opt_s){
	my $frac_v4_greater = $v4_has_more/$total;
	my $frac_unchanged  = $no_change/$total;
	my $frac_v4_lower   = $v4_has_less/$total;

	print "Fraction of V4 genes with more transcritps = $frac_v4_greater\n";
	print "Fraction of counts that are unchanged = $frac_unchanged\n";
	print "Fraction of V4 gene with fewer transcrpts = $frac_v4_lower\n";
    }
}
#-----------------------------------------------------------------------------
sub parse_gff{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	my @cols = split(/\t/, $line);
	if ($cols[2] eq 'mRNA' || 
	    $cols[2] eq 'transcript'){
	    if ($cols[8] =~ /Parent=gene:(\S+?);/|| 
		$cols[8] =~ /Parent=(\S+?);/){
		$DATA{$1}++;
	    }
	    else {die "you need to go back and fix the regex\n";}
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

