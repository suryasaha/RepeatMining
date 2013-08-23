#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Jun 25, 2012

use strict;
use warnings;
use Getopt::Long;
use IO::Handle;
use POSIX;
eval {
	require Bio::SeqIO;
};
use Bio::SeqIO;

=head1 NAME

preScreenREPCLASS.pl - Remove seqs with low complexity repeats from a repeat library file 

=head1 SYNOPSIS

  % preScreenREPCLASS.pl -l replib.fas -c 0.7 -m 100
  
=head1 DESCRIPTION

This script reads in a repeat library constructed by RepeatScout or another repeat finder and 
screens out consensus seqs with low complexity repeats > cutoff. Uses Tandem Repeat Finder with
params mostly from RepeatMasker. Use 5 instead of 7 for sensitive search. Writes out repeat rich
and short seqs in separate files.
          'matchWeight' => 2,
          'mismatchPenalty' => 5,
          'delta' => 5,
          'pm' => 80,
          'pi' => 10,
          'minScore' => 50,
          'maxPeriod' => 500,
          
=head1 NOTE
Make sure trf is in your path.

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --lib       <>  File with repeat library (required)
   --cutoff    <>  Max allowed % of low complexity repeats (float,suggested 0.70 fro 70%, required)
   --minlen    <>  Minimum length required (int,required)
      
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

my ($i,$j,$lib,$cutoff,$minlen,@temp,$ctr1,$ctr2);

GetOptions (
	'lib=s' => \$lib,
	'cutoff=f' => \$cutoff,
	'minlen=i' => \$minlen) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($lib) or (system('pod2text',$0), exit 1);
if (!(-e $lib)){print STDERR "$lib not found: $!\n"; exit 1;}
if($minlen < 0){
	system('pod2text',$0), exit 1;
}
if(($cutoff < 0) || ($cutoff >100)){
	system('pod2text',$0), exit 1;
}


my $in = Bio::SeqIO->new(-file=>$lib, -format=>'Fasta');
my $outFiltered = Bio::SeqIO->new(-file=>">filtered.${lib}", -format=>'Fasta');
my $outRep = Bio::SeqIO->new(-file=>">repeats.${lib}", -format=>'Fasta');
my $outShort = Bio::SeqIO->new(-file=>">short.${lib}", -format=>'Fasta');
$ctr1=$ctr2=0;
while (my $obj = $in->next_seq()){
	#chk length
	if($obj->length >= $minlen){
		my ($outTemp,@args,$rec,$ssrCumLen);
		#write out fasta with seq
		$outTemp = Bio::SeqIO->new(-file=>'>temp.fas', -format=>'Fasta');
		$outTemp->write_seq($obj);
		IO::Handle->autoflush();#explicitly flushing stream to disk before callin trf
		
		#run trf on it
		#trf NC_014774.fna 2 5 5 80 10 50 500 -h
		@args=("/home/surya/bin/trf",'temp.fas',2,5,5,80,10,50,500,'-h');
		#system(@args)==0 or die "system @args failed: $?";
		system(@args);#not checking for err, returns 256, ie 256/256 = 1 or err..dunno why
		
		#1374 1413 21 1.9 21 85 10 61 37 7 37 17 1.78 AGGAGATGAAGCAGGATTAAC AGGAGATGGAGCGGTATTAACAGGAGATGAAGCAGGATTA
		#6056 6103 22 2.2 22 78 10 58 41 6 8 43 1.60 TATTTAATTAAAAATATGATCA TATTTGATTAAAAATGTGATCATATTTAATTAAAACTATAGTCAATTT
		$i=`cat temp.fas.2.5.5.80.10.50.500.dat | wc -l`; chomp $i;
		$ssrCumLen=0;
		if($i>15){#ssr records are line 16+
			#read in coor ranges of repeats, compute sum
			unless(open(REP,"<temp.fas.2.5.5.80.10.50.500.dat")){print "not able to open temp.fas.2.5.5.80.10.50.500.dat\n\n";exit 1;}
			for(1..15){$rec=<REP>;}#skip first 16 lines
			while($rec=<REP>){#count length of all SSRs
				@temp=split(' ',$rec); $ssrCumLen+=$temp[1]-$temp[0]+1;
			}
		close(REP);
		}
		unlink('temp.fas.2.5.5.80.10.50.500.dat');
		
		#write to $outFiltered if below threshold
		if($ssrCumLen < ($cutoff*($obj->length))){ $outFiltered->write_seq($obj);}
		else{ $outRep->write_seq($obj); $ctr1++;}
		unlink("temp.fas");
	}
	else{
		$outShort->write_seq($obj); $ctr2++;
	}
}

if($ctr1==0){unlink("repeats.${lib}");}
if($ctr2==0){unlink("short.${lib}")}
exit;