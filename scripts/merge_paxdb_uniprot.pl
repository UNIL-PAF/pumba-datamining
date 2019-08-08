#! /usr/bin/perl

use strict;

####################
# merge paxdb abundance table with uniprot ids
#
# usage:
# ./merge_paxdb_uniprot.pl paxdb-uniprot-links-v4.1.tsv 9606/9606-HCT116_Koch_2014_emPAI.txt > paxdb_HCT116.txt


die "wrong parameters\n" unless scalar(@ARGV) == 2;


# print the header
print("protein_ac\tabundace\n");

open(LINKS, $ARGV[0]) or die "cannot open file [".$ARGV[0]."]\n";

my %links;

# parse links from paxdb ids to Uniprot ids
while(my $l = <LINKS>){
	chomp $l;
	my @a = split(/\s/, $l);
	$links{$a[0]} = $a[1];
}

close LINKS;

# parse abundace file and add corresponding Uniprot ids
open(ABUNDANCES, $ARGV[1]) or die "cannot open file [".$ARGV[1]."]\n";

while(my $l = <ABUNDANCES>){
	chomp $l;
	next if $l =~ /^#/;
	my @a = split(/\s/, $l);
	my $uniprot_ac = $links{$a[1]};
	
	# skip entries where we didnt find a match
	if(! $uniprot_ac){
		warn("could not find matching uniprot ac for [".$a[1]."]\n");
	 	next;
	}

	print($uniprot_ac."\t".$a[2]."\n");
}

close ABUNDANCES;



