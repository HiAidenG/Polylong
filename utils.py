## Collection of functions for my Pol II elongation project, refactoring all of these in biopython
from Bio import SeqIO
import pandas as pd
import numpy as np


def get_high_attention_kmers(AttMat, zScore = 1.0, positionCutoff = 100, freqCutoff = 0.75):
    """Get the highest scoring kmers from the attention matrix, return two lists of unique
    positive and negative kmers.
    
    === Args === 
    AttMat: pandas dataframe of the attention matrix. Rows are sequence kmers, columns are positions. 
    Values must be float. 
    
    === Filtering Parameters ===
    zScore: Selects kmers with z-scores >= zScore
    
    positionCutoff: Must be between 1 and final position in sequence inclusive.
    Selects kmers occuring in the first positionCutoff positions of the sequence. 
    
    freqCutoff: Selects kmers with frequencies in the top freqCutoff percentile.
    0 <= freqCutoff <= 1.0
    
    === Returns ===
    posKmers: list of unique positive kmers
    negKmers: list of unique negative kmers
    """
    