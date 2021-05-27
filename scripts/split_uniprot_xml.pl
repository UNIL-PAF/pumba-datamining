#! /usr/bin/perl
use strict;

# usage:
# ./split_uniprot_xml.pl /Users/rmylonas/tmp/datamining_pumba/uniprot/ < /Users/rmylonas/tmp/datamining_pumba/refprot_human.xml

my $out_path = $ARGV[0];

my $started = 0;
my $buffer;

while(my $l=<STDIN>){
	if(! $started){ 
		if($l =~ /\<entry/){
			$started = 1;
			$buffer  = $l;
		}
	}elsif($buffer && $l =~ /accession\>(.+)\</){
		open(OUT, ">", $out_path.$1.".xml") or die "Cannot open file in $out_path: $!\n";
		print OUT ($buffer.$l);
		$buffer = undef;
	}elsif($l =~ /\<\/entry/){
		print OUT ($l);
		close(OUT);
		$started = 0;
	}else{
		if(! $buffer){ 
			print OUT ($l);
		}else{
			$buffer .= $l;
		}
	}
}

