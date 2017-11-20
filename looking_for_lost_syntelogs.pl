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
my $LU = common_stuff::build_lu($FILE_1);

parse($FILE_2);

#PostData($LU);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my %data;

    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	
	my @cols = split(/\t/, $line);
	my $gid = 'mike';
	if ($cols[8] =~ /Parent=\S+?:(\S+?_F\S+?);/ || $cols[8] =~ /Parent=\S+?:(\S+?)_\S+?;/){
	    $gid = $1;
#	    print $gid,"\n";
	}
	my ($eid) = $cols[8] =~ /Name=(\S+?);/;
	my ($tid) = $cols[8] =~ /Parent=(\S+?);/;
	my $elen = $cols[4] - $cols[3] +1;
	my $overlap = $cols[-1];
	my $frac_ovl = $overlap/$elen;
	
	if (!defined($data{$gid}{$tid}{$eid}) || 
	    $data{$gid}{$tid}{$eid} < $frac_ovl){
	    $data{$gid}{$tid}{$eid} = $frac_ovl
	}
	else {next;}
	
	
    }
    $fh->close();
#    PostData(\%data);
    my %covered;
    foreach my $miss (keys %{$LU}){
	$miss =~ s/_\S+// unless $miss =~ /_F\S+/;
#	print $miss , "\n";
	foreach my $tid (keys %{$data{$miss}}){
#	print $miss, " $tid\n";

	    my $num_ex=0;
	    my $num_cov_ex=0;
	    foreach my $eid (keys %{$data{$miss}{$tid}}){
		$num_ex++;
		$num_cov_ex++ if $data{$miss}{$tid}{$eid} > .7;
	    }
	    my $frac_ex_covered = $num_cov_ex/$num_ex;
	    #if ($frac_ex_covered <= 0.7 && $frac_ex_covered > 0){
		if (defined($covered{$miss})){
		    if ($covered{$miss} > $frac_ex_covered){
			$covered{$miss} = $frac_ex_covered;
		    }
		    else{
			next;
		    }
		}
		else {
		    $covered{$miss} = $frac_ex_covered;
		}
		#print $miss ,"\n";
	    #}
	}
    }
    #PostData(\%covered);
    foreach my $gid (keys %covered){
	print "$gid\t$covered{$gid}\n";
    }
}
#-----------------------------------------------------------------------------

