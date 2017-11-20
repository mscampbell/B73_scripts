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

\n\n";

my $FILE_1 = $ARGV[0];
my $FILE_2 = $ARGV[1];
die($usage) unless $ARGV[0];

my $l_a = parse($FILE_1, $FILE_2);


#PostData($l_a);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file1 = shift;       
    my $file2 = shift;
    my $slu = common_stuff::build_lu($file1);
    my %aln_len;
    my %longest_aln;

    my $fh = new FileHandle;
    $fh->open($file2);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	my @cols = split(/\t/, $line);
	if ($cols[2] eq 'match_part'){
	    my ($pid) = $cols[8] =~ /Parent=(\S+?);/;
	    my ($target) = $cols[8] =~ /Target=(\S+?)\s/;
	    $target =~ s/_.*$//;
	    my ($len) = $cols[8] =~ /Length=(\S+?);/;
	    $aln_len{$target}{$pid} += $len;
	}
    }
    $fh->close();
    foreach my $hit (keys %aln_len){
	#print $hit ,"\n";
	foreach my $tlen (sort {$aln_len{$hit}{$b} <=> $aln_len{$hit}{$a}} keys %{$aln_len{$hit}}){
	    #I don't think I have to do anything here
	    $longest_aln{$hit}=$tlen unless defined($longest_aln{$hit});
	    #print "\t$tlen\t$aln_len{$hit}{$tlen}\n";
	}

    }
    
    foreach my $id (keys %longest_aln){
	#print "$id\n";
	print "$id\t$longest_aln{$id}\n" if defined($slu->{$id});
    }
    foreach my $x (keys %{$slu}){
	#print $x,"\n";
    }
    return (\%longest_aln);
}
#-----------------------------------------------------------------------------

