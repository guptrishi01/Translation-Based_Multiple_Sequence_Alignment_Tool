#!/usr/bin/env python3

"""
This module provides functions for sorting the sequences' ORFs 
based on k-mer similarity and create a guide tree for alignment

It includes functions for finding k-mers, calculating Jaccard Similarity for sorting,
and creating a guide tree based on the order. The goal is to use the constructed
guide tree for pairwise alignment.

Functions:
    - get_kmers(seq, k): Creates list of k-mers for sequence
    - build_binary_tree(order): Creates binary tree based on ordering of sequences
    - order_by_kmer_similarity(kmer_dict): Established order of sequences based on k-mer similarity
    - jaccard_calculation(kmers_1, kmers_2): Finds Jaccard similarity score of two list of k-mers
    - dataframe_to_dict(df, k): Creates dictionary of amino acid k-mers
"""

import pandas as pd
import sys

    # Remember which DNA fragment corresponds to LORF
    # After AA alignment, convert to nucleotide

def build_binary_tree(order: list) -> dict:
    """
    Construct a progressive binary tree from an ordered list of sequence labels.

    Args:
        order (list): List of sequences in order 

    Returns:
        terminals (dict): Dictionary mapping terminal node IDs to sequence labels
        internals (dict): Dictionary mappinginternal node IDs to children
    """
    # Create terminal nodes
    terminals = {i+1: name for i, name in enumerate(order)}
    n = len(order)

    # Error handling - one sequence means only one terminal node and no internal
    if n <= 1:
        return terminals, {}

    # Initialize internal nodes
    internal_nodes = {}
    next_id = n + 1

    # Create internal node of first two sequences
    internal_nodes[next_id] = [order[0], order[1]]
    prev = next_id
    next_id += 1

    # Create internal ndoes for rest of sequences
    for i in range(2, n):
        internal_nodes[next_id] = [prev, order[i]]
        prev = next_id
        next_id += 1

    return terminals, internal_nodes



def get_kmers(seq : str, k : int = 3) -> list:
    """
    Generate all k-mers of length k from a sequence.

    Args:
        seq (str): Amino acid sequence
        k (int) : Size of k-mer to create (default is 3)

    Returns:
        (list) : List of k-mers from sequence
    """

    if k < 1:
        sys.stderr.write("k must be a positive integer")
        raise ValueError("k must be a positive integer")
    

    return [seq[i: i + k] for i in range(len(seq) - k + 1)]

def jaccard_calculation(kmers_1 : list, kmers_2: list) -> float:
    """
    Compute the Jaccard similarity between two k-mer lists.

    Jaccard similarity = |intersection(k1, k2)| / |union(k1, k2)|.

    Args:
        kmers_1 (list): List of k_mers from first sequence
        kmers_2 (list): List of k_mers from second sequence

    Returns:
        (float) : Jaccard similarity value of 2 sequences
    """

    # No k-mers found - sequences are completely different
    if not kmers_1 or not kmers_2:
        return 0.0
    
    # Convert lists into sets and calculate intersection and union
    intersection = len(set(kmers_1) & set(kmers_2))
    union = len(set(kmers_1) | set(kmers_2))

    # Return 0 if either calculations are 0, else return jaccard value
    if union == 0 or intersection == 0:
        return 0.0
    
    return float(intersection / union)

def dataframe_to_dict(df : pd.DataFrame, k : int = 3) -> dict:
    """
    Convert a DataFrame of amino acid sequences into a dictionary of k-mers.

    Args:
        df (pd.DataFrame): Pandas Dataframe containing ORF amino acid sequences
        k (int): Length of k-mer

    Returns:
        kmer_dict (dict) : Dictionary mapping sequence IDs to list of k-mers for eachs sequence
    """

    kmer_dict = {}

    # Go through Pandas Dataframe and create dictionary entry {sequence id, list of k-mers}
    for idx, row in df.iterrows():
        amino_seq = str(row["Amino_Acids"])
        kmer_dict[str(row["id"])] = get_kmers(amino_seq, k)

    return kmer_dict

def order_by_kmer_similarity(kmer_dict: dict) -> list:
    """
    Produce a sequence ordering based on k-mer Jaccard similarity.

    Args:
        kmers_dict(dict) : Dictionary containing sequence IDs and list of k-mers for each sequence

    Returns:
        ordered (list) : List that contains order of sequences based on Jaccard similarity
    """
    # Obtain sequence IDs
    names = list(kmer_dict.keys())
    n = len(names)

    # If there is one sequence, return it
    if n <= 1:
        return names

    # Compute Jaccard similarities for each sequence pairing
    similarities = {}
    for i in range(n):
        for j in range(i+1, n):
            a, b = names[i], names[j]
            similarities[(a, b)] = jaccard_calculation(kmer_dict[a], kmer_dict[b])

    # Find the most similar pair
    best_pair = None
    best_val = -1.0
    for pair, v in similarities.items():
        if v > best_val:
            best_val = v
            best_pair = pair

    # Add most similar pair to list and remove it for consideration
    ordered = [best_pair[0], best_pair[1]]
    remaining = set(names) - set(ordered)

    # Iterate through rest of pairings and add them to list based on similarity
    while remaining:
        best_option = None
        best_score = -1.0

        for cand in remaining:
            # Grab pairing with highest similarity
            score = max(
                similarities.get((a, cand), similarities.get((cand, a), 0.0))
                for a in ordered
            )
            if score > best_score:
                best_score = score
                best_option = cand
        # Add pairing to list and remove it for consideration
        ordered.append(best_option)
        remaining.remove(best_option)

    return ordered
