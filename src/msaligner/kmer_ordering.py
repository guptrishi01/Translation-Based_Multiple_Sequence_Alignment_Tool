#!/usr/bin/env python3

"""
This program generates a distance matrix of amino acid ORF sequences sorted by k-mer similarity

Author: Rishi Gupta
Contact: rgupta25@charlotte.edu

"""

import pandas as pd
import math
from itertools import combinations


    # Remember which DNA fragment corresponds to LORF
    # After AA alignment, convert to nucleotide

def build_binary_tree(order: list):
    """
    Build a progressive binary tree from ordering.
    Returns:
        terminals : {1: label, 2: label, ...}
        internals : {nodeId: [leftChild, rightChild]}
    """
    terminals = {i+1: name for i, name in enumerate(order)}
    n = len(order)

    if n <= 1:
        return terminals, {}

    internal_nodes = {}
    next_id = n + 1

    # First join
    internal_nodes[next_id] = [order[0], order[1]]
    prev = next_id
    next_id += 1

    # Add remaining sequences
    for i in range(2, n):
        internal_nodes[next_id] = [prev, order[i]]
        prev = next_id
        next_id += 1

    return terminals, internal_nodes



def get_kmers(seq : str, k : int = 4) -> list:
    return [seq[i: i + k] for i in range(len(seq) - k + 1)]

def jaccard_calculation(kmers_1 : list, kmers_2: list) -> float:
    """
    Calculate the Jaccard similarity score for two sequences
    Args:
        seq1 (str): The first amino acid sequence
        seq2 (str): The second amino acid sequence
        length (int): The length of k-mers specified by user (default is 4)

    Returns:
        float: Mash Distance of two amino acid sequences
    """

    # No k-mers found - sequences are completely different
    if not kmers_1 or not kmers_2:
        return 0.0
    

    intersection = len(set(kmers_1) & set(kmers_2))
    union = len(set(kmers_1) | set(kmers_2))

    if union == 0 or intersection == 0:
        return 0.0
    
    return float(intersection / union)

def dataframe_to_dict(df : pd.DataFrame, k = 4) -> dict:
    kmer_dict = {}

    for idx, row in df.iterrows():
        amino_seq = str(row["Amino_Acids"])
        kmer_dict[str(row["id"])] = get_kmers(amino_seq, k)

    return kmer_dict

def order_by_kmer_similarity(kmer_dict: dict) -> list:
    names = list(kmer_dict.keys())
    n = len(names)

    if n <= 1:
        return names

    # Precompute pairwise similarities
    sim = {}
    for i in range(n):
        for j in range(i+1, n):
            a, b = names[i], names[j]
            sim[(a, b)] = jaccard_calculation(kmer_dict[a], kmer_dict[b])

    # 1) Find most similar pair
    best_pair = None
    best_val = -1.0
    for pair, v in sim.items():
        if v > best_val:
            best_val = v
            best_pair = pair

    ordered = [best_pair[0], best_pair[1]]
    remaining = set(names) - set(ordered)

    # 2) Add sequences most similar to the cluster
    while remaining:
        best_option = None
        best_score = -1.0

        for cand in remaining:
            score = max(
                sim.get((a, cand), sim.get((cand, a), 0.0))
                for a in ordered
            )
            if score > best_score:
                best_score = score
                best_option = cand

        ordered.append(best_option)
        remaining.remove(best_option)

    return ordered
