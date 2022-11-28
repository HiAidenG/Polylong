import numpy as np
from scipy.stats import zscore

class kMerDF:
    """
    This is a special class designed for working with transformer
    attention scores, specifically with DNA sequences that have
    been parsed into k-mers. Designed to work alongside the 
    kMer class.

    === Instance Variables ===
    @np.array df: a numpy array of kMer objects

    === Representation Invariants ===
    
    === Preconditions ===
    - scores is a numpy array of attention scores
    - dnadf is a numpy array of DNA sequences
    - labels is a numpy array of labels

    scores.shape == dnadf.shape
    """

    def __init__(self, scores, dnadf, labels):
        scores = zscore(scores)
        self.df = self._make_kmerdf(scores, dnadf, labels)
    
    def _build_kmer(self, ID, pos, score, seq, label):
        """
        Initializes a kMer object and returns it.
        """
        return kMer(ID, pos, score, seq, label)

    def _make_kmerdf(self, scores, dnadf, labels):
        """
        Creates a numpy dataframe of kMer objects from the scores and
        dnadf dataframes. The kMer objects are stored in a pandas
        dataframe for easy access and manipulation.
        """
        kmerdf = []
        for row_idx, row_val in enumerate(scores):
            row = []
            for col_idx, _ in enumerate(row_val):
                ID = row_idx
                pos = col_idx
                score = round(scores[row_idx][col_idx], 3)
                seq = dnadf[row_idx][col_idx]
                label = labels[row_idx]
                kmer = self._build_kmer(ID, pos, score, seq, label)
                row.append(kmer) 
            kmerdf.append(row)
        return np.array(kmerdf)

    def __len__(self):
        """
        Returns the length of the kMerDF object.
        """
        return len(self.df)

    def __getitem__(self, key):
        """
        Returns the kMer object at the given key.
        """
        return self.df[key]
            
    def __iter__(self):
        """
        Returns an iterator over the kMerDF object.
        """
        return kMerDF_iter(self)

    def head(self, n):
        """
        Returns the first n rows of the kMerDF object.

        Returns the empty list if n = 0. Returns all rows of the kMerDF
        if n > len(self.df).
        """
        return self.df[:n] 
        
    def filter_kmers(self, zscore=0, seq_freq=0, pos_cutoff=0):
        """
        Filters the k-mer dataframe and returns a list of kMer objects
        satisfying the given parameters.

        === Preconditions ===
        - float zscore: the z-score threshold, default 0
        - int seq_freq: the minimum frequency of a sequence, default 0
                        0 < seq_freq < 100
        - int pos_cutoff: the maximum position of the k-mer, after which subsequent
                            k-mers are ignored, default 0
                          0 < pos_cutoff < len(self.df[0])
        """
        if pos_cutoff > 0:
            r_df = self.df[: , :pos_cutoff]
        if zscore > 0:
            r_df = r_df.filter_by_zscore(zscore)
        if seq_freq > 0:
            r_df = r_df.filter_by_seq_freq(seq_freq)
        return r_df

    def filter_by_zscore(self, zscore):
        """
        Mutates a kMerDF object by removing all kMers below the given zscore threshold.
        """
        for row in self.df:
            for kmer in row:
                if kmer.score < zscore:
                    row.remove(kmer)
        return self

    def filter_by_seq_freq(self, seq_freq_percentile):
        """
        Mutates a kMerDF object by removing all kMers below the given sequence frequency threshold.

        === Precondition ===
        - int seq_freq: the minimum frequency of a sequence, default 0
                        0 < seq_freq_percentile < 100
        """
        for row in self.df:
            for kmer in row:
                if self.get_sequence_frequency(kmer.seq) < seq_freq_percentile:
                    row.remove(kmer)
        return self

        # TODO: figure out how to get percentiles for sequence frequency


    def get_sequence_frequency(self, seq):
        """
        Returns the frequency of a given sequence in the kMerDF object.
        """
        pass
       
        

    

    def find_matching_seqs(self, seq):
        """
        Returns a list of kMer objects the match the given sequence

        === Preconditions ===
        - str seq: the sequence to search for
        """
        kmers = []
        for row in self.df:
            for kmer in row:
                if kmer.seq == seq:
                    kmers.append(kmer)
        return kmers

    # TODO: this is buggy, it returns a list of empty lists, which I don't want
    # Probably want to split this up into some smaller methods. 
    def get_contiguous_kmers(self, zscore):
        """
        Returns List[List[kMer]], where kMers satisfy the zscore 
        threshold and each List[kMer] is a contiguous high-scoring sequence.
        """
        return_list = []
        for i in range(len(self.df)):
            contiguous_kmers = []
            row = self.df[i]
            for j in range(1, len(row)-1): #start at index one, end at second to last index (to avoid index out of bounds)
                kmer = row[j]
                if kmer.score >= zscore and (row[j+1].score >= zscore or row[j-1].score >= zscore): # if the next kmer is also high scoring
                    contiguous_kmers.append(kmer) # add it to the list
            # edge cases, check the first and last kmers
            if row[0].score >= zscore and row[1].score >= zscore:
                contiguous_kmers.append(row[0])
            if row[-1].score >= zscore and row[-2].score >= zscore:
                contiguous_kmers.append(row[-1])
            return_list.append(contiguous_kmers) 
        return return_list

class kMerDF_iter:
    """
    Iterator for kMerDF class.
    """
    def __init__(self, kmerdf):
        self.kmerdf = kmerdf
        self.row = 0
        self.col = 0
        
    def __next__(self):
        if self.row >= len(self.kmerdf):
            raise StopIteration
        else:
            kmer = self.kmerdf[self.row][self.col]
            self.col += 1
            if self.col >= len(self.kmerdf[self.row]):
                self.col = 0
                self.row += 1
            return kmer

class kMer:
    """
    Designed for working with transformer attention matrices.
    Works in conjunction with the attMat class.

    === Instance Variables ===
    @str ID: the identifer for the DNA sequence this k-mer belongs to. 
    @int pos: the position of this k-mer in the sequence
    @int score: the attention score of this k-mer
    @string seq: the sequence of this k-mer.
    @int label: the classficiation label for this k-mer. 0 for unstable, 1 for stable.

    === Representation Invariants ===
    - 0 <= pos <= len(seq) - k
    - k > 0
    - seq is a valid DNA sequence
    """

    def __init__(self, ID, pos, score, seq, label):
        self.ID = ID
        self.pos = pos
        self.score = score
        self.seq = seq
        self.label = label

    def __str__(self):
        return f"Sequence ID: {self.ID}, Sequence: {self.seq}, Position: {self.pos}, Score: {self.score}"  

#function for retrieving the unique sequences from a list of kMers
def get_unique_seqs(kmers):
    """
    Returns a list of unique sequences from a list of kMers. 
    Chooses a representative kMer, the first to appear in the list. 
    
    You can search for all kmers matching the representative's sequence
    using kMerDF.search(<kMer.seq>) 
    """
    seqs = []
    for kmer in kmers:
        match = False
        for seq in seqs: 
            if kmer.seq == seq.seq:
                match = True
        if not match:
            seqs.append(kmer)
    return seqs

def kmer2seq(kmers):
    """
    Convert kmers to original sequence
    
    Arguments:
    kmers -- str, kmers separated by space.
    
    Returns:
    seq -- str, original sequence.
    """
    kmers_list = kmers.split(" ")
    bases = [kmer[0] for kmer in kmers_list[0:-1]]
    bases.append(kmers_list[-1])
    seq = "".join(bases)
    assert len(seq) == len(kmers_list) + len(kmers_list[0]) - 1
    return seq

    
