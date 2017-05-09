from Bio import SeqIO
import sys
from multiprocessing import Pool

def generate_wildcard_kmers(kmer):
	wc_kmers = []
	
	for i in xrange(len(kmer)):
		new_kmer = kmer[0:i] + "N" + kmer[i+1:]
		wc_kmers.append(new_kmer)
	
	for i in xrange(len(kmer)-1):
		new_kmer = kmer[0:i] + "NN" + kmer[i+2:]
		wc_kmers.append(new_kmer)
	
	return wc_kmers
	
def predict_label_given_seq(seq):
	score = 0
	k = 8
	for i in xrange(len(seq)-k+1):
		kmer = seq[i:i+k]
		if kmer in weights:
			score+=weights[kmer]
		
		for wc_kmer in generate_wildcard_kmers(kmer):
			if wc_kmer in weights:
				score+=weights[wc_kmer]

	return 1 if score>0 else -1
	
def predict_label_and_print(record):
	mut_seq_label = predict_label_given_seq(record.seq)
	print record.id+"\t"+str(mut_seq_label)

if __name__=="__main__":
	w_file = sys.argv[1]
	orig_seqs_file = sys.argv[2]
	mut_seqs_file = sys.argv[3]
	
	weights = dict()
	with open(w_file, 'r') as f:
		for line in f:
			splits = line.strip().split()
			kmer = splits[0]
			weight = float(splits[1])
			weights[kmer] = weight
		
	orig_seq_labels = dict()
	orig_seq_label_file = open("seql_orig_seq_labels.txt", 'w')
	with open(orig_seqs_file, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			orig_seq_label = predict_label_given_seq(record.seq)
			orig_seq_label_file.write(record.id+"\t"+str(orig_seq_label)+"\n")
	orig_seq_label_file.close()		
	
	pool = Pool(32)
	records = []
	with open(mut_seqs_file, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			records.append(record)
	
	pool.map(predict_label_and_print, records)