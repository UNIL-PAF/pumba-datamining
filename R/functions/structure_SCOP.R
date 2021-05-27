# parse the SCOP database information

# load database
csv_file_path <- "/Users/rmylonas/Work/PAF/projects/pumba/data/datamining/scop/190917_scop.csv"
scop_classes_path <- "/Users/rmylonas/Work/PAF/projects/pumba/data/datamining/scop/SCOP_classes.csv"

scop_table <- read.csv(csv_file_path, sep="\t")
scop_classes <- read.csv(scop_classes_path)

# get scop classes for a given AC
get_scop_classes <- function(protein_ac){
  as.character(scop_table$scopClasses[scop_table$ac == protein_ac])
}