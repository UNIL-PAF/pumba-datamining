# used to get the charge for each protein
fasta_seq_path <- "/Users/admin/Work/PAF/projects/pumba/data/fasta/20170609_UP000005640_9606.csv"
fasta_seqs <- read.csv(file=fasta_seq_path)

# check if a sequence exists
seq_exists <- function(protein_ac){ sum(fasta_seqs$ac == protein_ac) == 1 }


# get the sequence
get_seq <- function(protein_ac){ as.character(fasta_seqs$seq[fasta_seqs$ac == protein_ac]) }
