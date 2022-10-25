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
    

# test the attMat class
# make sure init works

scores = [[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
          [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]]

def test_attMat_init()
    