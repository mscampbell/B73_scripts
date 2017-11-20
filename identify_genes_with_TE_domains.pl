#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;
use gff3_annotation_stuff;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
this script takes the interproscan results from the ensembl pipeline and matches them up with the IPR ids form josh.
\n\n";

my $FILE_FROM_JOSH   = $ARGV[0]; #contianed IDS 
my $FILE_FROM_ANDREW = $ARGV[1]; #contains pfam ids for transcripts
my $FILE_FROM_SHARON = $ARGV[2]; 
my $GFF3_FILE        = $ARGV[3];
die($usage) unless $ARGV[2];
my %JOSH;
my %ANDREW;
my %SHARON;
my %IPR2PFAM;
#my ($GFF3_HR, $GID_TID_LU_HR) = gff3_annotation_stuff::parse_gff($GFF3_FILE);
parse($FILE_FROM_JOSH, 'j');
parse($FILE_FROM_ANDREW, 'a');

#parse($FILE_FROM_SHARON, 's');
report_te_pfams();

#PostData(\%IPR2PFAM);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report_te_pfams{
    foreach my $x (keys %IPR2PFAM){
	print join("\n", @{$IPR2PFAM{$x}}),"\n" if defined($JOSH{$x});
	
    }
#    PostData(\%JOSH);
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my $flag = shift;
    my $fh = new FileHandle;
    my %printed;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	if ($flag eq 'j'){
	   my @x = split(/\s/, $line);
	   $JOSH{$x[0]}=1;
	} 
	elsif ($flag eq 'a'){
	   next if $line =~ /^id/;
	   my @x = split(/\t/, $line);
	   my $ipr = shift(@x);
	   foreach my $y (@x){
	       if ($y =~ /^\S/){
		   #print "$ipr\n";
		   #print $y , "\n";
		   push(@{$ANDREW{$y}}, $ipr);
		   push(@{$IPR2PFAM{$ipr}}, $y);
		    
	       }
	   }
	} 
	elsif ($flag eq 's'){
	    my @x = split(/\t/, $line);
	    my $tid = $x[0];
	    my $did = $x[4];
	    foreach my $converted_id (@{$ANDREW{$did}}){
		if (defined($JOSH{$converted_id})){
		    print "$tid\n" unless defined($printed{$tid});
		    $printed{$tid} = 1;
		    
                    #print "$tid\t$converted_id\t$did\n";
		    #push(@{$SHARON{$tid}{$converted_id}}, $did);
		}
	    }
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

