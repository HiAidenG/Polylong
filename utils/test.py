# Simple test suite for the utils module

import kmers
from kmers import kMer, kMerDF

# Test kMer class
def test_kmer():
    kmer = kMer(id = 1, seq = "ATG", score = 0.5, pos = 1, label=1)
    assert kmer.id == 1
    assert kmer.seq == "ATG"
    assert kmer.score == 0.5
    assert kmer.pos == 1
    assert kmer.label == 1

# Test kMerDF class
def test_kmerdf():
    
    