from Bio import SeqIO

#load in sequence data
FA = "../data/raw_data.fa"
SEQ_LENGTH = 700

#clean up the sequence data
def clean_seqs(FA, SEQ_LENGTH):
    sequences = {}
    for record in SeqIO.parse(FA, "fasta"):
        sequence = str(record.seq).upper()
        #remove mitochondrial sequences and sequences of improper length
        if 'chrM' not in record.id and len(record.seq) == SEQ_LENGTH:
            if record.id not in sequences:
                sequences[record.id] = sequence
    
    with open('../data/clean_seqs.fa', 'w') as f:
        for key in sequences:
            f.write(">" + key + "\n" + sequences[key] + "\n")


#split data into training and testing sets
#using chromosomes 1-16 for training and rest for testing
def split_seqs(clean_fasta):
    training = ['chr' + str(i) for i in range(1, 17)]
    for record in SeqIO.parse(clean_fasta, "fasta"):
        for chromosome in training:
            if chromosome in record.id:
                with(open('../data/training_seqs.fa', 'a')) as f:
                    f.write(">" + str(record.id) + "\n" + str(record.seq) + "\n")
            else:
                with(open('../data/testing_seqs.fa', 'a')) as f:
                    f.write(">" + str(record.id) + "\n" + str(record.seq) + "\n")

if __name__ == "__main__":
    clean_seqs(FA, SEQ_LENGTH)
    split_seqs("../data/clean_seqs.fa")
