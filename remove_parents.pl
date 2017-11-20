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
my %LU;
parse($FILE, 'lu');
parse($FILE, 'ex');


#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my $tag = shift;

    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	next if $line =~ /^#/;
	my @col = split(/\t/, $line);
	if ($tag eq 'lu'){
	    if ($col[2] eq 'mRNA'){
		my ($id) = $col[8] =~ /ID=(\S+?);/;
		$LU{$id}=1;
	    }
	}
	if ($tag eq 'ex'){
	    if ($col[2] eq 'exon'){
		my ($p_list) = $col[8] =~ /Parent=(.*)/;
		my @pars = split(/\,/, $p_list);
		foreach my $par (@pars){
		    if (defined($LU{$par})){
			$line =~ s/Parent=\S+/Parent=$par/;
			
		    }
		}
	    }
	    print $line."\n";
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------
    
