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
parse($FILE, 'rp');
#parse($FILE, 'rn');
#PostData(\%LU);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my $flag = shift;
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	my @cols = split(/\t/, $line);
	if ($flag eq 'lu'){
	    if ($cols[2] eq 'mRNA'){
		#removing the -RA remember to put it back on for the
		#non-gene featurs
		#New development the -RA might not be -RA
		my ($gname) = $cols[8] =~ /Name=(\S+?)-\S+?;/;
		my ($mname) = $cols[8] =~ /Name=(\S+?);/;
		my ($gid) = $cols[8] =~ /Parent=(\S+?);/;
		my ($tid) = $cols[8] =~ /ID=(\S+?);/;
		$LU{'g'}{$gid} = $gname;
		$LU{'m'}{$tid} = $mname
	    }
	}
	#elsif ($flag eq 'rn'){
	#    #added this as a hack to fix some names in step 3 annotaiton of B73
	#    
	#    if($cols[2] eq 'mRNA'){
	#	my ($tid) = $cols[8] =~ /ID=(\S+?);/;
	#	$line =~ s/Parent=\S+?;/Parent=$LU{'g'}{$gid};/;
	#	$line =~ s/ID=\S+?;/ID=$LU{'m'}{$tid};/;
	#	print $line , "\n";
	#    }
	#}
	elsif ($flag eq 'rp'){
	    if ($cols[2] eq 'gene'){
		my ($id) = $cols[8] =~ /ID=(\S+?);/;
		$line =~ s/$id/$LU{'g'}{$id}/g;
		print $line , "\n";
	    }	
	    elsif($cols[2] eq 'mRNA'){
		my ($gid) = $cols[8] =~ /Parent=(\S+?);/;
		my ($tid) = $cols[8] =~ /ID=(\S+?);/;
		$line =~ s/Parent=\S+?;/Parent=$LU{'g'}{$gid};/;
		$line =~ s/ID=\S+?;/ID=$LU{'m'}{$tid};/;
		print $line , "\n";
	    }
	    elsif($cols[2] eq 'exon' ||
		  $cols[2] eq 'five_prime_UTR' ||
		  $cols[2] eq 'three_prime_UTR' ||
		  $cols[2] eq 'CDS'){
		#fix the exon ID
		my @ats = split(/;/, $cols[8]);
		my ($id) = $ats[0] =~ /ID=(\S+?):/;
		$ats[0] =~ s/ID=$id/ID=$LU{'m'}{$id}/;
		my $new_col9 = $ats[0].";Parent=";
		#now work on the parents
		$ats[1] =~ s/Parent=//;
		my @pts = split(/\,/, $ats[1]);
		foreach my $pt (@pts){
		    $new_col9 .= $LU{'m'}{$pt} . ",";
		}
		$new_col9 =~ s/,$//;
		$cols[8] = $new_col9;
		
		print join("\t", @cols), "\n";
		
		
	    }
	    else {print STDERR "WARNING: you missed this feature \"$cols[2]\"\n"}
	}    
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

