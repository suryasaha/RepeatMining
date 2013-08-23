#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Feb 5, 2011

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use Switch;

=head1 NAME

trfReportSummary.pl - Generate summary staatistics for Tandem Repeat Finder report 

=head1 SYNOPSIS

  % trfReportSummary.pl  -r report.dat
  
=head1 DESCRIPTION

TRF command trf seqs.fa 2 5 5 80 10 50 500 -d -h 1>/dev/null
Generate following stats
Total number of identified exact SSRs (inexact not implemented) 
Number of sequences containing exact di, tri, tetra, penta, all SSRs
Number of sequences containing inexact di, tri, tetra, penta, all SSRs
Number of sequences containing more than 1 exact/inexact SSR  (not implemented)
Number of sequences without SSRs
Total number of seqs

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --report    <>  Report (required)
      
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

my ($i,$j,$rep,@temp,$ctr,%stats,$rec,$flagA);

GetOptions (
	'report=s' => \$rep) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($rep) or (system('pod2text',$0), exit 1);
if (!(-e $rep)){print STDERR "$rep not found: $!\n"; exit 1;}

$stats{'seq'}=0; $stats{'ssr'}=0; $stats{'sssr'}=0; $stats{'dissr'}=0; $stats{'disssr'}=0; $stats{'trissr'}=0; 
$stats{'trisssr'}=0; $stats{'tetssr'}=0; $stats{'tetsssr'}=0; $stats{'restssr'}=0; $stats{'restsssr'}=0; 
$stats{'seqnossr'}=0;

unless(open(REP,"<$rep")){print "not able to open $rep\n\n";exit 1;}
for(1..8){$rec=<REP>;}#skip first 8 lines
$flagA=0;
while($rec=<REP>){
	#seq name
	if($rec=~ /^Sequence\: /){
		$stats{'seq'}++; 
		if($flagA==1){ $stats{'seqnossr'}++;}
		$flagA=1;
	}
	
	#ssr
	if($rec=~ /^[1-9]/){
		$flagA=0; #reset
		$stats{'ssr'}++;
		
		#process records
		@temp=split(" ",$rec);
		switch($temp[2]){
			case 2 {
				if($temp[5]==100){$stats{'dissr'}++;}
				else {$stats{'disssr'}++;}
			}
			case 3 {
				if($temp[5]==100){$stats{'trissr'}++;}
				else {$stats{'trisssr'}++;}
			}
			case 4 {
				if($temp[5]==100){$stats{'tetssr'}++;}
				else {$stats{'tetsssr'}++;}
			}
			else {
				if($temp[5]==100){$stats{'restssr'}++;}
				else {$stats{'restsssr'}++;}
			}
		}
	}
}
close(REP);

#print summary of report
print "Sequence summary\nSequences\t",$stats{'seq'},"\n";
print "Sequences with SSRs\t",$stats{'seq'}-$stats{'seqnossr'},"\n";
print "Sequences without SSRs\t",$stats{'seqnossr'},"\n\nSSR summary\n";
print "Dinucleotide exact SSRs\t",$stats{'dissr'},"\n";
print "Dinucleotide inexact SSRs\t",$stats{'disssr'},"\n";
print "Trinucleotide exact SSRs\t",$stats{'trissr'},"\n";
print "Trinucleotide inexact SSRs\t",$stats{'trisssr'},"\n";
print "Tetranucleotide exact SSRs\t",$stats{'tetssr'},"\n";
print "Tertanucleotide inexact SSRs\t",$stats{'tetsssr'},"\n";
print "Rest exact SSRs\t",$stats{'restssr'},"\n";
print "Rest inexact SSRs\t",$stats{'restsssr'},"\n";

exit;