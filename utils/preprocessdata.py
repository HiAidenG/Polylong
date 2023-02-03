from Bio import SeqIO
import random
import argparse
import os

#Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--pos_file', type=str, help='Path to positive sequences')
parser.add_argument('-n', '--neg_file', type=str, help='Path to negative sequences')
parser.add_argument('-t', '--type', type=bool, default=False, help='Whether the sequences should be processed into windows for training multiple models, default is False. NOTE: This assumes the sequences are 1024 bp long, behavior is undefined if they are not. This is also a bit of a file dump, consider redirecting output.')
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
            print("Expected length: " + str(length), "Actual length: " + str(len(record.seq)))
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

def extract_window(sequences: dict, window_start: int, window_end: int) -> dict:
    """
    Returns a subset of the given dictionary, containing the nucleotides between the window_start and window_end indices.
    """
    subset = {}
    for ID, seq in sequences.items():
        subset[ID] = (seq[0][window_start:window_end], seq[1])
    return subset

def kmerize_seqs(sequences: dict) -> None:
    """
    Mutates the given dictionary by converting the sequences to kmers
    """
    for ID, seq in sequences.items():
        sequences[ID] = (seq2kmer(seq[0], 6), seq[1])
              

def shuffle_seqs(sequences: dict) -> None:
    """
    Shuffle the positions of the sequences in the dictionary, such that the classes aren't written in order
    Returns None because the dictionary is mutated
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
        f.write('sequence\tlabel\n')
        for ID, seq in sequences.items():
            f.write(seq[0] + '\t' + str(seq[1]) + '\n')

def single_window() -> tuple:
    positive_records = list(SeqIO.parse(POS_FILE, 'fasta'))
    negative_records = list(SeqIO.parse(NEG_FILE, 'fasta'))

    #Clean sequences
    positives = clean_seqs(positive_records)
    negatives = clean_seqs(negative_records)

    #Reformat negative sequences
    #negatives = reformat_seqs(negatives)

    #Merge positive and negative sequences
    merged = merge_seqs(positives, negatives)

    #kmerize sequences
    kmerize_seqs(merged)

    #Shuffle sequences
    shuffle_seqs(merged)

    #Split sequences into training, validation, and testing sets
    train, dev, test = split_seqs(merged)
    
    return train, dev, test

def tiling_pipeline(window_start: int, window_end: int) -> tuple:
    # There's a lot of duplicate code here, but I think it makes it more explicit what's going on
    positive_records = list(SeqIO.parse(POS_FILE, 'fasta'))
    negative_records = list(SeqIO.parse(NEG_FILE, 'fasta'))

    #Clean sequences
    positives = clean_seqs(positive_records)
    negatives = clean_seqs(negative_records)

    #Reformat negative sequences
    negatives = reformat_seqs(negatives)

    #Extract window
    positives = extract_window(positives, window_start, window_end)
    negatives = extract_window(negatives, window_start, window_end)

    #Merge positive and negative sequences
    merged = merge_seqs(positives, negatives)

    #kmerize sequences
    kmerize_seqs(merged)

    #Shuffle sequences
    shuffle_seqs(merged)

    #Split sequences into training, validation, and testing sets
    train, dev, test = split_seqs(merged)
    
    return train, dev, test


def tiled_window():
    # Again, lots of duplicate code here... this one I could probably refactor a bit
    # One full length sequence of 1024 bp:
    # 0-1024
    (train, dev, test) = single_window()
    os.mkdir('0-1024')
    write_seqs(train, '0-1024/train.tsv')
    write_seqs(dev, '0-1024/dev.tsv')
    write_seqs(test, '0-1024/test.tsv')

    # Three overlapping windows of 512 bp:
    # 0-512, 256-768, 512-1024
    (train, dev, test) = tiling_pipeline(0, 512)
    os.mkdir('0-512')
    write_seqs(train, '0-512/train.tsv')
    write_seqs(dev, '0-512/dev.tsv')
    write_seqs(test, '0-512/test.tsv')

    (train, dev, test) = tiling_pipeline(256, 768)
    os.mkdir('256-768')
    write_seqs(train, '256-768/train.tsv')
    write_seqs(dev, '256-768/dev.tsv')
    write_seqs(test, '256-768/test.tsv')

    (train, dev, test) = tiling_pipeline(512, 1024)
    os.mkdir('512-1024')
    write_seqs(train, '512-1024/train.tsv')
    write_seqs(dev, '512-1024/dev.tsv')
    write_seqs(test, '512-1024/test.tsv')

    # Seven overlapping windows of 256 bp:
    # 0-256, 128-384, 256-512, 384-640, 512-768, 640-896, 768-1024
    (train, dev, test) = tiling_pipeline(0, 256)
    os.mkdir('0-256')
    write_seqs(train, '0-256/train.tsv')
    write_seqs(dev, '0-256/dev.tsv')
    write_seqs(test, '0-256/test.tsv')

    (train, dev, test) = tiling_pipeline(128, 384)
    os.mkdir('128-384')
    write_seqs(train, '128-384/train.tsv')
    write_seqs(dev, '128-384/dev.tsv')
    write_seqs(test, '128-384/test.tsv')

    (train, dev, test) = tiling_pipeline(256, 512)
    os.mkdir('256-512')
    write_seqs(train, '256-512/train.tsv')
    write_seqs(dev, '256-512/dev.tsv')
    write_seqs(test, '256-512/test.tsv')

    (train, dev, test) = tiling_pipeline(384, 640)
    os.mkdir('384-640')
    write_seqs(train, '384-640/train.tsv')
    write_seqs(dev, '384-640/dev.tsv')
    write_seqs(test, '384-640/test.tsv')

    (train, dev, test) = tiling_pipeline(512, 768)
    os.mkdir('512-768')
    write_seqs(train, '512-768/train.tsv')
    write_seqs(dev, '512-768/dev.tsv')
    write_seqs(test, '512-768/test.tsv')

    (train, dev, test) = tiling_pipeline(640, 896)
    os.mkdir('640-896')
    write_seqs(train, '640-896/train.tsv')
    write_seqs(dev, '640-896/dev.tsv')
    write_seqs(test, '640-896/test.tsv')

    (train, dev, test) = tiling_pipeline(768, 1024)
    os.mkdir('768-1024')
    write_seqs(train, '768-1024/train.tsv')
    write_seqs(dev, '768-1024/dev.tsv')
    write_seqs(test, '768-1024/test.tsv')

if __name__ == '__main__':
    if args.type == False:
        (train, dev, test) = single_window()
        write_seqs(train, 'train.tsv')
        write_seqs(dev, 'dev.tsv')
        write_seqs(test, 'test.tsv')
    else:
        tiled_window()