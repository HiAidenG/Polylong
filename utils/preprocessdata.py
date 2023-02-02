from Bio import SeqIO
import random
import argparse 

#Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--pos_file', type=str, help='Path to positive sequences')
parser.add_argument('-n', '--neg_file', type=str, help='Path to negative sequences')
parser.add_argument('-t', '--type', type=bool, default=False, help='Whether the sequences should be processed into windows for training multiple models, default is False')
args = parser.parse_args()

POS_FILE = args.pos_file
NEG_FILE = args.neg_file

""" 
This script exists to preprocess data to be used for training the testing a DNABERT model. 

DNABERT requires tsv files with the following format:
    - Each line is a sequence with a label
    - Sequence and label are tab separated
    - Label is 0 or 1
    - Sequence is a string of kmers separated by spaces
    
By default, this script will take in two fasta files, one for the positive sequences and one for the negative sequences.
It will clean the sequences by removing any mitochondrial genes and ensuring all sequences are of the same length.

It will then split the sequences into training and testing sets, using chromosomes 1-14 for training, 15-16 for validation, and the rest for testing.
Sequences will be converted into kmers and written to their respective tsv files. The script outputs test.tsv, train.tsv, and dev.tsv to the data directory.

TODO: type flag to split sequences into windows for training multiple models
"""

#Update: 1/31/23: I'm trying to reduce I/O operations, so I'm changing the functions to do everything in memory

def clean_seqs(records: list) -> dict:
    """
    Return a dictionary of cleaned sequences from a list of SeqIO records.
    Cleaned sequences are all the same length and do not contain mitochondrial genes.

    Issues a warning if there is a sequence length mismatch.
    """
    sequences = {}
    length = len(records[0].seq)
    for record in records:
        sequence = str(record.seq).upper()
        if 'chrM' not in record.id and len(record.seq) == length:
            if record.id not in sequences:
                sequences[record.id] = sequence
        elif len(record.seq) != length:
            print('Warning: sequence length mismatch')
    return sequences

def seq2kmer(seq, k) -> str:
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

def reformat_seqs(sequences: dict) -> dict:
    """
    Negatives sequences are formatted as chr:start-end, so this function reformats them to be chr_start-end
    returns a dictionary of reformatted sequences
    """
    reformatted = {}
    for ID, seq in sequences.items():
        temp = ID.split(':')
        ID = str(temp[2]) + '_' + str(temp[3])
        reformatted[ID] = seq
    return reformatted

def kmerize(sequences: dict) -> dict:
    """
    Convert sequences to kmers
    """
    kmers = {}
    for ID, seq in sequences.items():
        kmers[ID] = seq2kmer(seq, 6)
    return kmers

def merge_seqs(positives: dict, negatives: dict) -> dict:
    """
    Merge positive and negative sequences into a single dictionary, keeping track of the labels

    Returns a dictionary where each key is the sequence ID and each value is a tuple of the sequence and the label
    """
    merged = {}
    for ID, seq in positives.items():
        merged[ID] = (seq, 1)
    for ID, seq in negatives.items():
        merged[ID] = (seq, 0)
    return merged

def kmerize_seqs(sequences: dict) -> None:
    """
    Mutates the given dictionary by converting the sequences to kmers
    """
    for ID, seq in sequences.items():
        label = seq[1]
        seq = seq2kmer(seq[0], 6)
        sequences[ID] = (seq, label)
    
def shuffle_seqs(sequences: dict) -> None:
    """
    Mutates the given dictionary by shuffling the order of the sequences
    """
    keys = list(sequences.keys())
    random.shuffle(keys)
    shuffled = {}
    for key in keys:
        shuffled[key] = sequences[key]
    sequences = shuffled

def split_seqs(sequences: dict) -> tuple:
    """
    Returns a tuple of dictionaries. Split into training, validation, and testing sets

    Training set is chromosomes 1-14
    Validation set is chromosomes 15-16
    Testing set is chromosomes 17-22
    """
    train = {}
    dev = {}
    test = {}
    for ID, seq in sequences.items():
        #yes, this is very ugly
        if 'chr1' in ID or 'chr2' in ID or 'chr3' in ID or 'chr4' in ID or 'chr5' in ID or 'chr6' in ID or 'chr7' in ID or 'chr8' in ID or 'chr9' in ID or 'chr10' in ID or 'chr11' in ID or 'chr12' in ID or 'chr13' in ID or 'chr14' in ID:
            train[ID] = seq
        elif 'chr15' in ID or 'chr16' in ID:
            dev[ID] = seq
        else:
            test[ID] = seq
    return train, dev, test

def write_seqs(sequences: dict, write_dir: str):
    """
    Write sequences to a tsv file

    At this point, I don't care about sequence IDs, so I'm just writing the sequence and the label
    """
    with open(write_dir, 'w') as f:
        for ID, seq in sequences.items():
            f.write(seq[0] + '\t' + str(seq[1]) + '\n')

def main():
    positive_records = list(SeqIO.parse(POS_FILE, 'fasta'))
    negative_records = list(SeqIO.parse(NEG_FILE, 'fasta'))


    #Clean sequences
    positives = clean_seqs(positive_records)
    negatives = clean_seqs(negative_records)

    #Reformat negative sequences
    # negatives = reformat_seqs(negatives)

    #Merge positive and negative sequences
    merged = merge_seqs(positives, negatives)

    print(merged)

    #kmerize sequences
    merged = kmerize_seqs(merged)

    print(merged)


    #Split sequences into training, validation, and testing sets
    # train, dev, test = split_seqs(merged)


    #Write sequences to tsv files
    # write_seqs(train, 'train.tsv')
    # write_seqs(dev, 'dev.tsv')
    # write_seqs(test, 'test.tsv')


if __name__ == '__main__':
    main()

    



























# #additional function for reformatting sequences
# def reformat_seqs(fasta):
#     with open('neg_reformat.fa', 'w') as f:
#         for record in SeqIO.parse(fasta, "fasta"):
#             temp = record.id.split(':')
#             record.id = str(temp[2]) + '_' + str(temp[3])
#             f.write(">" + record.id + "\n" + str(record.seq) + "\n")



# def get_kmers(pos_fa, neg_fa, to_write):
#     #Read in fasta file
#     positives = []
#     negatives = []
#     for seq_record in SeqIO.parse(pos_fa, 'fasta'):
#         if len(seq_record) == SEQ_LENGTH:
#             seq = str(seq_record.seq)
#             positives.append(seq)

#     for seq_record in SeqIO.parse(neg_fa, 'fasta'):
#         if len(seq_record) == SEQ_LENGTH:
#             seq = str(seq_record.seq)
#             negatives.append(seq)

#     with open(to_write, 'w') as f:
#         for record in positives:
#             record = seq2kmer(record, 6)
#             f.write("%s\t%s\n" % (record,'1'))
#         for record in negatives:
#             record = seq2kmer(record, 6)
#             f.write("%s\t%s\n" % (record,'0'))

#     # Shuffle the sequences
#     with open(to_write, 'r') as f:
#         lines = f.readlines()
#         random.shuffle(lines)
#     with open(to_write, 'w') as f2:
#         f2.write("%s\t%s\n" % ('sequence','label'))
#         f2.writelines(lines)
    

# if __name__ == "__main__":
#     import pandas as pd
#     clean_seqs(POS_FILE, SEQ_LENGTH)
#     split_seqs('pos_clean_seqs.fa')

#     clean_seqs(NEG_FILE, SEQ_LENGTH)
#     #The UU sequences have to be reformatted to match the format of the positive sequences.
#     reformat_seqs('neg_clean_seqs.fa')
#     split_seqs('neg_reformat.fa')

#     get_kmers('pos_training_seqs.fa', 'neg_training_seqs.fa', 'train.tsv')

#     get_kmers('pos_testing_seqs.fa', 'neg_testing_seqs.fa', 'test.tsv')

#     df = pd.read_csv('test.tsv', sep='\t')
#     dev_set = df[-1000:]
#     df = df[:-1000]
#     df.to_csv('test.tsv', sep='\t', index=False)
#     dev_set.to_csv('dev.tsv', sep='\t', index=False)

#     print("Done!")
