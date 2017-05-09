from Bio import SeqIO
import numpy as np
import sys
import argparse


def generate_wildcard_kmers(kmer, max_consec_wildcard):
    if 0<max_consec_wildcard and max_consec_wildcard>len(kmer):
        sys.exit("max_consec_wildcard not in range 0 to k")
    wc_kmers = []
    
    for num_consec in xrange(1,max_consec_wildcard+1):
        for i in xrange(len(kmer)-num_consec+1):
            new_kmer = kmer[0:i] + "N"*num_consec + kmer[i+num_consec:]
            wc_kmers.append(new_kmer)
            
    return wc_kmers


def kmerize_fa(input_fasta, k, max_consec_wildcard):
    kmer_vocab=set()
                 
    for kmer in kmer_vocab:
        if not len(kmer) == k:
            sys.exit("Input kmer vocab should contain kmers of length k")
    
    all_kmer_counts = dict()
    
    numSeqs = 0
    for record in SeqIO.parse(input_fasta, 'fasta'):
        kmer_counts = kmerize_seq(record.seq, k, max_consec_wildcard)
        all_kmer_counts[record.id] = kmer_counts
        kmer_vocab.update(kmer_counts.keys())
        numSeqs+=1
    
    
    kmer_vocab = list(kmer_vocab)
    design_matrix = np.zeros((numSeqs, len(kmer_vocab)), dtype=np.int)
    
    i=0
    for record_id in all_kmer_counts:
        kmer_counts = all_kmer_counts[record_id]
        for (j,kmer) in enumerate(kmer_vocab):
            if kmer in kmer_counts:
                count = kmer_counts[kmer]
                design_matrix[i][j] = count            
        i+=1 
    
    return design_matrix, kmer_vocab
    
def kmerize_seq(seq, k, max_consec_wildcard):
    if 0<max_consec_wildcard and max_consec_wildcard>k:
        sys.exit("max_consec_wildcard not in range 0 to k")
                 
    kmer_counts = dict()
    for i in xrange(len(seq)-k+1):
        kmer = str(seq[i:i+k])
        if kmer not in kmer_counts:
            kmer_counts[kmer] = 0
        
        kmer_counts[kmer] = kmer_counts[kmer] + 1
        
        for wc_kmer in generate_wildcard_kmers(kmer, max_consec_wildcard):
            if wc_kmer not in kmer_counts:
                kmer_counts[wc_kmer] = 0
            kmer_counts[wc_kmer] = kmer_counts[wc_kmer] + 1

    return kmer_counts
                 
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Kmerizer for Fasta file')
    parser.add_argument('-i', '--fasta', help='fasta file input', required=True)
    parser.add_argument('-k', '--kmer-length', type=int, help='length of k-mer desired', required=True)
    parser.add_argument('-m', '--max-consec-wc', type=int, help='maximum continuous wildcards', required=True)
    
    args = parser.parse_args()
    design_matrix, kmer_vocab = kmerize_fa(args.fasta, args.kmer_length, args.max_consec_wc)
    
    np.savetxt('design_matrix.txt',design_matrix, delimiter = "\t")
    
    with open('kmer_vocab.txt', 'w') as vocab_out:
        for word in kmer_vocab:
            vocab_out.write(word+"\n")
    