from Bio import SeqIO
import random

#TODO: this script needs a lot of work, I should probably rewrite the whole thing

#load in sequence data
POS_FILE = '../data/posd512.fa'
NEG_FILE = '../data/negd512.fa'
SEQ_LENGTH = 512

#clean up the sequence data
def clean_seqs(fasta, SEQ_LENGTH):
    sequences = {}
    for record in SeqIO.parse(fasta, "fasta"):
        sequence = str(record.seq).upper()
        #remove mitochondrial sequences and sequences of improper length
        if 'chrM' not in record.id and len(record.seq) == SEQ_LENGTH:
            if record.id not in sequences:
                sequences[record.id] = sequence
    if fasta == POS_FILE:
        write_to = 'pos_clean_seqs.fa'
    elif fasta == NEG_FILE:
        write_to = 'neg_clean_seqs.fa'
    else:
        print('Error: invalid fasta file')
    with open(write_to, 'w') as f:
        for key in sequences:
            f.write(">" + key + "\n" + sequences[key] + "\n")

#split data into training and testing sets
#using chromosomes 1-16 for training and rest for testing
def split_seqs(clean_fasta):
    if clean_fasta == 'pos_clean_seqs.fa':
        write_to = 'pos_training_seqs.fa'
        write_to2 = 'pos_testing_seqs.fa'
    elif clean_fasta == 'neg_reformat.fa':
        write_to = 'neg_training_seqs.fa'
        write_to2 = 'neg_testing_seqs.fa'
    else:
        print('Error: invalid fasta file')
    # training = ['chr' + str(i) for i in range(1, 17)]
    for record in SeqIO.parse(clean_fasta, "fasta"):
        #find the chromosome number for record
        chr_num = record.id.split('_')[0]
        if chr_num == 'chr1' or chr_num == 'chr2' or chr_num == 'chr3' or chr_num == 'chr4' or chr_num == 'chr5' or chr_num == 'chr6' or chr_num == 'chr7' or chr_num == 'chr8' or chr_num == 'chr9' or chr_num == 'chr10' or chr_num == 'chr11' or chr_num == 'chr12' or chr_num == 'chr13' or chr_num == 'chr14' or chr_num == 'chr15' or chr_num == 'chr16':
            #if the chromosome number is in the preset training set, write to training file
                with open(write_to, 'a') as f:
                    f.write(">" + record.id + "\n" + str(record.seq) + "\n")
        # if the chromosome number is not in the preset training set, write to testing file
        else:
            with open(write_to2, 'a') as f:
                f.write(">" + record.id + "\n" + str(record.seq) + "\n")

#additional function for reformatting sequences
def reformat_seqs(fasta):
    with open('neg_reformat.fa', 'w') as f:
        for record in SeqIO.parse(fasta, "fasta"):
            temp = record.id.split(':')
            record.id = str(temp[2]) + '_' + str(temp[3])
            f.write(">" + record.id + "\n" + str(record.seq) + "\n")

def seq2kmer(seq, k):
    #function from DNABERT, link to paper: https://academic.oup.com/bioinformatics/article-abstract/37/15/2112/6128680
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

def get_kmers(pos_fa, neg_fa, to_write):
    #Read in fasta file
    positives = []
    negatives = []
    for seq_record in SeqIO.parse(pos_fa, 'fasta'):
        if len(seq_record) == SEQ_LENGTH:
            seq = str(seq_record.seq)
            positives.append(seq)

    for seq_record in SeqIO.parse(neg_fa, 'fasta'):
        if len(seq_record) == SEQ_LENGTH:
            seq = str(seq_record.seq)
            negatives.append(seq)

    with open(to_write, 'w') as f:
        for record in positives:
            record = seq2kmer(record, 6)
            f.write("%s\t%s\n" % (record,'1'))
        for record in negatives:
            record = seq2kmer(record, 6)
            f.write("%s\t%s\n" % (record,'0'))

    # Shuffle the sequences
    with open(to_write, 'r') as f:
        lines = f.readlines()
        random.shuffle(lines)
    with open(to_write, 'w') as f2:
        f2.write("%s\t%s\n" % ('sequence','label'))
        f2.writelines(lines)
    

if __name__ == "__main__":
    import pandas as pd
    clean_seqs(POS_FILE, SEQ_LENGTH)
    split_seqs('pos_clean_seqs.fa')

    clean_seqs(NEG_FILE, SEQ_LENGTH)
    #The UU sequences have to be reformatted to match the format of the positive sequences.
    reformat_seqs('neg_clean_seqs.fa')
    split_seqs('neg_reformat.fa')

    get_kmers('pos_training_seqs.fa', 'neg_training_seqs.fa', 'train.tsv')

    get_kmers('pos_testing_seqs.fa', 'neg_testing_seqs.fa', 'test.tsv')

    df = pd.read_csv('test.tsv', sep='\t')
    dev_set = df[-1000:]
    df = df[:-1000]
    df.to_csv('test.tsv', sep='\t', index=False)
    dev_set.to_csv('dev.tsv', sep='\t', index=False)

    print("Done!")
