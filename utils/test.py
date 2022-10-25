from attMat import kMer, attMat

# Module for testing the attMat and kMer classes.

def test_kmer():
    # Test the kMer class
    kmer = kMer(1, 10, 0.5, "ATCG", 0)
    assert kmer.id == 1
    assert kmer.pos == 10
    assert kmer.score == 0.5
    assert kmer.seq == "ATCG"
    assert kmer.label == 0
    
    
    