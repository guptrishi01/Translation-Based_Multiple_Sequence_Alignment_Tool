#!/usr/bin/env python3

"""
This program generates a distance matrix of amino acid ORF sequences sorted by k-mer similarity

Author: Rishi Gupta
Contact: rgupta25@charlotte.edu

"""

import pandas as pd
import numpy as np
from Bio import SeqIO
import math
from itertools import combinations

def kmer_similarity(df: pd.DataFrame, length=4) -> list:
    """
    Create a distance matrix to order amino acid sequences by k-mer similarity
    Args:
        df (pd.DataFrame): Pandas DataFrame containing amino acid ORF sequences
        length (int): The length of k-mers specified by user (default is 4)

    Returns:
        list: List containing order of sequences
    """

    """ Dictionary for Amino Acids """
    kmer_array = {}
    
    """ Lists that will contain kmers for similarity calculations """

    # Generate all unique two-digit combinations
    two_digit_combinations = list(combinations(range(len(df['Amino_Acids'])), 2))
    # Go through each combination and calculate Mash Distance of sequences
    for combination in two_digit_combinations:
        kmer_array[(combination[0],combination[1])] = index_and_distance(df['Amino_Acids'][combination[0]],df['Amino_Acids'][combination[1]] ,length)
    
    # Sort the dictionary by ascending Mash Distance values
    keys, values = dict(sorted(kmer_array.items(), key=lambda item: item[1]))
    # Return pairs of sequences to be aligned based on similarity (most to least)
    return keys



# ATGAGAGAGAGA
## [ATG, TGA, GAG, ..., AGA]
## Do this for a pair of sequences 
## For kmer in list of kmers in seq 1 if kmer in list of kmers in seq2 count ++
## Find the common k-mers (How many k-mers they have in common - an integer like 55) - dictates similarity
### More kmers they share - more similar they are
### Do that for all your sequences and keep track of scores
## Whatevers your most similar gets aligned first, and you follow the scores

def index_and_distance(seq1 : str, seq2: str, length: int) -> float:
    kmers_1 = [seq1[i: i + length] for i in range(len(seq1) - length + 1)]
    kmers_2 = [seq2[i: i + length] for i in range(len(seq2) - length + 1)]
    jaccard = len(set(kmers_1) & set(kmers_2)) / len(set(kmers_1) | set(kmers_2))
    # No k-mers found - sequences are completely different
    if jaccard == 0:
        return 1
    return float((-1 / length) * math.log((2 * jaccard) / (1 + jaccard)))