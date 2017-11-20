#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use recover_IR;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
I moved the recover_IR.pl script into a module to increase flexibility and abstract out the problems that I already solved. This driver passes clusters of gene annotations and isoseq data to the recover_IR.pm module as an array reference of the lines from a gff3 file. 
\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];
#someware around here I need to collect the clusters before I send them to the parse subroutien. I am going to see if tabix indexing can be used to make it fast, but I will probably start with generating an array of arrays just to make sure it works on a first pass and judge performance then.
my $gene_regions = get_gene_regions($FILE, 200);
my $clusters= get_clusters($gene_regions, $FILE); # gets an array ref of array refs
run_recover_IR($clusters);

#my $file_array = parse($FILE);
#recover_IR::run($file_array);
#PostData($clusters);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub run_recover_IR{
    print STDERR "\nsending the clusters to the recover_IR module\n\n";
    my $clusters = shift;
    my $count = 0;
    my $num_clusters = @$clusters;
    while (defined(@{$clusters->[0]})){
	my $cluster = shift @{$clusters};
	#PostData($cluster);
	recover_IR::run($cluster);
	$count++;
	if ($count/1000 == int($count/1000)){
	    print STDERR "\nprocessed $count of $num_clusters clusters\n\n";
	    }
    }
}
#-----------------------------------------------------------------------------
sub get_clusters{
    print STDERR "\ngetting gene clusters ie the isoseq data in a gene region\n\n";
    my $gene_regions = shift;
    my $file = shift;
    my $num_gr = @$gene_regions;
    my @clusters;
    while (defined(@{$gene_regions->[0]})){
	my @region_cluster;
	my $region = shift @{$gene_regions};
	my ($sid, $beg, $end, $str) = @{$region};
	my $num_left = @$gene_regions;
	if($num_left/100 == int($num_left/100)){
	    print STDERR "\nprocessed", $num_gr - $num_left, " of $num_gr gene regions\n\n";
	}
	my $fh = new FileHandle;
	$fh->open($file); # this is probably slowing things down a lot
	my $start_time = time;
	print STDERR "\nstarted getting a cluster\n\n";
	while (defined(my $line = <$fh>)){
	    chomp($line);
	    last if $line =~ /\#\#FASTA/;
	    next if $line =~ /\#/;
	    
	    my @col = split(/\t/, $line);
	    if ($col[0] eq $sid && $col[6] eq $str && 
		($col[1] eq 'maker' || $col[1] eq 'est2genome:iso') &&
		($col[3] >= $beg && $col[4] <= $end)){
		push @region_cluster, $line;
	    }
	    else {next;}
#	push @file_array, $line;
	}
	push @clusters, \@region_cluster;
	my $end_time = time;
	print STDERR "\nfinished getting a cluster ", $end_time - $start_time, " seconds\n\n";

	$fh->close();
	
    }
        print STDERR "\nfinished getting gene clusters\n\n";
    return (\@clusters);
}
#-----------------------------------------------------------------------------
sub get_gene_regions{
    print STDERR "\ngetting gene regions\n\n";
    my $file = shift;       
    my $flank = shift;
    my @clusters;
    my @gene_regions;

    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /\#\#FASTA/;
	next if $line =~ /\#/;
	
	my @col = split(/\t/, $line);
	my $count=1;
	if ($col[2] eq 'gene'){
	    push (@gene_regions, [$col[0], ($col[3] - $flank), ($col[4] + $flank), $col[6]]);
	    $count++;
	}
#	push @file_array, $line;
    }
    print STDERR "\nfinished getting gene regions\n\n";
    $fh->close();
#    PostData(\@gene_regions);
    my $num_gene_regions = @gene_regions;
    print STDERR "\nidentified $num_gene_regions gene regions\n\n";
    return(\@gene_regions);
}
#-----------------------------------------------------------------------------
#sub parse{
#
#    my $file = shift;       
#    my @file_array;
#    my $fh = new FileHandle;
#    $fh->open($file);
#    
#    while (defined(my $line = <$fh>)){
#	chomp($line);
#	push @file_array, $line;
#    }
#    $fh->close();
#    return(\@file_array);
#}
#-----------------------------------------------------------------------------

