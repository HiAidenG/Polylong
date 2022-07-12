from Bio import SeqIO
import random

SEQ_LENGTH = 700
POS_FILE = '../data/pos_testing_seqs.fa'
NEG_FILE = '../data/neg_testing_seqs.fa'

TO_WRITE = 'dev.tsv'

# Read in the fasta file
test_positives = []
test_negatives =[]
for seq_record in SeqIO.parse(POS_FILE, 'fasta'):
    if len(seq_record) == SEQ_LENGTH:
        seq = str(seq_record.seq)
        test_positives.append(seq)

for seq_record in SeqIO.parse(NEG_FILE, 'fasta'):
    if len(seq_record) == SEQ_LENGTH:
        seq = str(seq_record.seq)
        test_negatives.append(seq)

def seq2kmer(seq, k):
    """
    Convert original sequence to kmers
    
    Arguments:
    seq -- str, original sequence.
    k -- int, kmer of length k specified.
    
    Returns:
    kmers -- str, kmers separated by space

    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

with open(TO_WRITE, 'w') as f:
    for record in test_positives:
        record = seq2kmer(record, 6)
        f.write("%s\t%s\n" % (record,'1'))
    for record in test_negatives:
        record = seq2kmer(record, 6)
        f.write("%s\t%s\n" % (record,'0'))


# Shuffle the sequences
with open(TO_WRITE, 'r') as f:
    lines = f.readlines()
    random.shuffle(lines)
    with open(TO_WRITE, 'w') as f2:
        f2.write("%s\t%s\n" % ('sequence','label'))
        f2.writelines(lines)

    






