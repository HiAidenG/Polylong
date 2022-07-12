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
    
    with open('../data/clean_pos_seqs.fa', 'w') as f:
        for key in sequences:
            f.write(">" + key + "\n" + sequences[key] + "\n")

#split data into training and testing sets
#using chromosomes 1-16 for training and rest for testing
def split_seqs(clean_fasta):
    # training = ['chr' + str(i) for i in range(1, 17)]
    for record in SeqIO.parse(clean_fasta, "fasta"):
        #find the chromosome number for record
        chr_num = record.id.split('_')[0]
        if chr_num == 'chr1' or chr_num == 'chr2' or chr_num == 'chr3' or chr_num == 'chr4' or chr_num == 'chr5' or chr_num == 'chr6' or chr_num == 'chr7' or chr_num == 'chr8' or chr_num == 'chr9' or chr_num == 'chr10' or chr_num == 'chr11' or chr_num == 'chr12' or chr_num == 'chr13' or chr_num == 'chr14' or chr_num == 'chr15' or chr_num == 'chr16':
            #if the chromosome number is in the preset training set, write to training file
                with open('../data/training_seqs.fa', 'a') as f:
                    f.write(">" + record.id + "\n" + str(record.seq) + "\n")
        # if the chromosome number is not in the preset training set, write to testing file
        else:
            with open('../data/testing_seqs.fa', 'a') as f:
                f.write(">" + record.id + "\n" + str(record.seq) + "\n")

#additional function for reformatting sequences
def reformat_seqs(FA):
    for record in SeqIO.parse(FA, 'fasta'):
        temp = record.id.split(':')
        record.id = str(temp[2]) + '_' + str(temp[3])
        with open('../data/reformat.fa', 'a') as f:
                    f.write(">" + record.id + "\n" + str(record.seq) + "\n")

if __name__ == "__main__":
    clean_seqs(FA, SEQ_LENGTH)
    split_seqs("../data/clean_pos_seqs.fa")
