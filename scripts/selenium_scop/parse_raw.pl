#! /usr/bin/perl


my $res_path = "/Users/admin/Work/PAF/projects/pumba/data/datamining/scop/orig/";
my $csv_ac_file = "/Users/admin/Work/PAF/projects/pumba/data/fasta/20170609_UP000005640_9606.csv";
my $threshold = 0.001;

open (CSV, $csv_ac_file) or die "canooot open file\n";

my %scop_classes = (	
	"a" => "All alpha proteins",
	"b" => "All beta proteins",
	"c" => "Alpha and beta proteins (a/b)",
	"d" => "Alpha and beta proteins (a+b)",
	"e" => "Multi-domain proteins (alpha and beta)",
	"f" => "Membrane and cell surface proteins and peptides",
	"g" => "Small proteins",
	"h" => "Coiled coil proteins",
	"i" => "Low resolution protein structures",
	"j" => "Peptides",
	"k" => "Designed proteins",
	"l" => "Artifacts"
);


# print the header of the result
print "ac\tscopClasses\n";

# remove header
<CSV>;


while(<CSV>){
	my %sel_classes;

	if(/(^.+?),.+/){
		my $ac = $1;

		# try to open and parse the corresponding raw file from webscrape_scop.py
		my $res_file = $res_path.$ac.".txt";
		if(-f $res_file){
			open (RES, $res_file) or die "Cannot open res file";	
			my $write_res = 0;

			while(my $l = <RES>){
				chomp $l;
				if($l =~ /Sequences producing significant alignments/){
					$write_res = 1;
				}elsif($l =~ /^>.+/){
					$write_res = 0;
				}elsif($write_res){
					if($l =~ /^\w+\s+(\w).+\s+([\w|\-|\.]+)/){
						$sel_classes{$1} = 1 if($2 <= $threshold);
					}
				}

			}

			close RES;

			print("$ac\t".join(',', sort(keys %sel_classes))."\n");

		}else{
			print STDERR "could not find file for [$ac]\n";
		}
		

	}
	
}



close CSV;
