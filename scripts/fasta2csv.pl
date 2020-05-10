#! /usr/bin/perl -w

use strict;

## USAGE
# ./fasta2csv.pl < /Users/admin/Work/PAF/projects/pumba/data/fasta/20170609_UP000005640_9606.fasta > /Users/admin/Work/PAF/projects/pumba/data/fasta/20170609_UP000005640_9606.csv


# change line seperator
$/="\n>";

# print the header
print("ac,seq,header\n");

while(my $l = <STDIN>){
	my @a = split("\n", $l);

	my $ac;
	if($a[0] =~ /^.+?\|(.+?)\|.+/){
		$ac = $1;
	}

	my $str = "";

	foreach my $i (1..$#a){
		last if($a[$i] =~ />/);
		$str .= $a[$i];
	}

	my $header = $a[0];
	$header =~ s/,/;/g;

	print("$ac,$str,".$header."\n") if($ac);	
}
