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
#PostData($LU);
parse($FILE_2);

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       

    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	my @cols = split(/\t/, $line);
	
	if($cols[1] =~ /^augustus/ || 
	   $cols[1] =~ /^fgene/ || 
	   $cols[1] =~ /^pred_gff/){
	    my $id ='mike';
	    if($cols[2] eq 'match'){
		($id) = $cols[8] =~ /Name=(\S+?);/;
		#print $id,"\n";
	    }
	    elsif($cols[2] eq 'match_part'){
		($id) = $cols[8] =~ /Target=(\S+)/;
	    }
	    print "$line\n" if defined($LU->{$id});
	    
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

