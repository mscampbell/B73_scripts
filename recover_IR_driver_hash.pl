#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use recover_IR;
use recover_IR_all;
use Getopt::Std;
use vars qw($opt_r $opt_a $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('ragpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
I moved the recover_IR.pl script into a module to increase flexibility and abstract out the problems that I already solved. This driver passes clusters of gene annotations and isoseq data to the recover_IR.pm module as an array reference of the lines from a gff3 file.
-a calls one module that gets all the changes
-r calls another module that only gets retentions of already annotated gene models
\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];
#someware around here I need to collect the clusters before I send them to the parse subroutien. I am going to see if tabix indexing can be used to make it fast, but I will probably start with generating an array of arrays just to make sure it works on a first pass and judge performance then.
my ($gene_regions, $num_gr) = get_gene_regions($FILE, 200); 
# struct looks like this
# $gene_regions{'gene_regions'}{scaffold_id} = arrary of arrays
# $gene_regions{'isoseq'}      {scaffold_id} = arrary of gff3 lines

my $clusters= get_clusters($gene_regions, $FILE, $num_gr); # gets an array ref of array refs
run_recover_IR($clusters);

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
	recover_IR::run($cluster) if $opt_r; #original that only captures retentions
	recover_IR_all::run($cluster) if $opt_a; #new version that captures all changes
	$count++;
	    print STDERR "\nprocessed $count of $num_clusters clusters\n\n";
    }
}
#-----------------------------------------------------------------------------
sub get_clusters{
    print STDERR "\ngetting evidence clusters ie the isoseq data in a gene region\n\n";
    my $gene_regions = shift;
    my $file = shift;
    #PostData($gene_regions);
# struct looks like this
# $gene_regions{'gene_regions'}{scaffold_id} = arrary of arrays
# $gene_regions{'isoseq'}      {scaffold_id} = arrary of gff3 lines
    my $num_gr = shift;
    my @clusters;
    my $count = 0;
    foreach my $scaff (keys %{$gene_regions->{'gene_regions'}}){
	while (defined($gene_regions->{'gene_regions'}->{$scaff}->[0])){
	    my @region_cluster;
	    my $region = shift @{$gene_regions->{'gene_regions'}->{$scaff}};
	    my ($sid, $beg, $end, $str) = @{$region};
	    my $num_left = $num_gr - $count;
	    if($num_left/100 == int($num_left/100)){ 
		print STDERR "\nprocessed ", $num_gr - $num_left, " of $num_gr gene regions\n\n";
	    }
	    
	    foreach my $line (@{$gene_regions->{'isoseq'}->{$scaff}}){
		#go through the isoseq and assign to a cluster
		my @col = split(/\t/, $line);
		if ($col[6] eq $str && ($col[3] >= $beg && $col[4] <= $end)){
		    push @region_cluster, $line;
		}
	    }
	    foreach my $line (@{$gene_regions->{'maker'}->{$scaff}}){
		#go through the genes and assign to a cluster
		my @col = split(/\t/, $line);
		if ($col[6] eq $str && ($col[3] >= $beg && $col[4] <= $end)){
		    push @region_cluster, $line;
		}
		
		else {next;}
	    }
	    $count++;
	    push @clusters, \@region_cluster;
	} 
    }
        print STDERR "\nfinished getting evidence clusters\n\n";
    return (\@clusters);
}
#-----------------------------------------------------------------------------
sub get_gene_regions{
    print STDERR "\ngetting gene regions\n\n";
    my $file = shift;       
    my $flank = shift;
    my @clusters;
    #my @gene_regions;
    my %gene_regions;
    my $num_gene_regions=0;
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /\#\#FASTA/;
	next if $line =~ /\#/;
	
	my @col = split(/\t/, $line);
	my $count=1;
	if ($col[2] eq 'gene'){
	    $num_gene_regions++;
	    push (@{$gene_regions{'gene_regions'}{$col[0]}}, [$col[0], ($col[3] - $flank), ($col[4] + $flank), $col[6]]);
	    $count++;
	}
	if ($col[1] eq 'est2genome:iso'){
	    push (@{$gene_regions{'isoseq'}{$col[0]}}, $line);
	}
	if ($col[1] eq 'maker'){
	    push (@{$gene_regions{'maker'}{$col[0]}}, $line);
	}
    }
    print STDERR "\nfinished getting gene regions\n\n";
    $fh->close();

    print STDERR "\nidentified $num_gene_regions gene regions\n\n";
    return(\%gene_regions, $num_gene_regions);
}
#-----------------------------------------------------------------------------
