# used to get the charge for each protein
fasta_seq_path <- "/Users/rmylonas/Work/PAF/projects/pumba/data/fasta/UP000005640_canonical_plus_splices_16092020_PUMBA.csv"
fasta_seqs <- read.csv(file=fasta_seq_path)

# check if a sequence exists
seq_exists <- function(protein_ac){ sum(fasta_seqs$ac == protein_ac) == 1 }


# get the sequence
get_seq <- function(protein_ac){ as.character(fasta_seqs$seq[fasta_seqs$ac == protein_ac]) }

# get theorethical mass
get_theo_mass <- function(protein_ac, database){
  seq_get <- getURL(paste0("localhost:9000/sequence/", protein_ac, "/database/", database))
  my_seq <- RJSONIO::fromJSON(seq_get)
  my_seq$molWeight
}


