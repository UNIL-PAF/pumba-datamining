#! /usr/bin/perl -w

use strict;

## USAGE
# ./fasta2csv.pl < /Users/rmylonas/Work/PAF/projects/pumba/data/fasta/20170609_UP000005640_9606.fasta > /Users/admin/Work/PAF/projects/pumba/data/fasta/20170609_UP000005640_9606.csv


# change line seperator
$/="\n>";

# print the header
print("ac,seq,header,protein_name,gene_name\n");

while(my $l = <STDIN>){
	my @a = split("\n", $l);

	my $ac;
	my $gene_name = '';
	my $protein_name;

	if($a[0] =~ /^.+?\|(.+?)\|.+?\s(.+?)\sOS.+GN=(.+?)\s.*/){
		$ac = $1;
		$protein_name = $2;
		$gene_name = $3;
	}elsif($a[0] =~ /^.+?\|(.+?)\|.+?\s(.+?)\sOS.*/){
		$ac = $1;
		$protein_name = $2;
	}

	my $str = "";

	foreach my $i (1..$#a){
		last if($a[$i] =~ />/);
		$str .= $a[$i];
	}

	my $header = $a[0];
	$header =~ s/,/;/g;
	if($protein_name){ $protein_name =~ s/,/;/g; }

	print("$ac,$str,".$header.",$protein_name,$gene_name\n") if($ac);	
}
