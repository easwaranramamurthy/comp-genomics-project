import numpy as np
from Bio import SeqIO

if __name__=="__main__":
	with open('yy2_motif_matches_ids.txt') as f:
		ids = set([line.strip() for line in f])
		
	seq_subset = []
	with open("ENCFF330LVI_seq.fa", "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			if record.id in ids:
				seq_subset.append(record)
	
	SeqIO.write(seq_subset, "ENCFF330LVI_matching_YY1_seq.fa", "fasta")
