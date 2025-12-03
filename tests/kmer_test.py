import pytest
import pandas as pd
from msaligner.kmer_ordering import (
    get_kmers,
    jaccard_calculation,
    dataframe_to_dict,
    order_by_kmer_similarity,
    build_binary_tree
)

def test_get_kmers_basic():
    """
    Testing correctness of get_kmers method
    """
    assert get_kmers("ABCDE", 3) == ["ABC", "BCD", "CDE"]

def test_get_kmers_short_sequence():
    """
    Testing invalid case of get_kmers method
    """
    assert get_kmers("AB", 3) == []

def test_get_kmers_invalid_k():
    """
    Testing check of ValueError
    """
    with pytest.raises(ValueError):
        get_kmers("ABC", 0)


def test_jaccard_basic():
    """
    Testing correctness of jaccard similarity calculation
    """
    k1 = ["AAA", "BBB", "CCC"]
    k2 = ["BBB", "CCC", "DDD"]
    assert jaccard_calculation(k1, k2) == 2/4

def test_jaccard_no_overlap():
    """
    Testing invalid case of jaccard similarity calculation
    """
    assert jaccard_calculation(["AAA"], ["BBB"]) == 0.0

def test_dataframe_to_dict_basic():
    """
    Testing correct creation of kmer dictionary
    """
    df = pd.DataFrame({
        "id": ["S1", "S2"],
        "Amino_Acids": ["ABCDE", "AAA"]
    })
    d = dataframe_to_dict(df, k=2)
    assert d["S1"] == ["AB", "BC", "CD", "DE"]
    assert d["S2"] == ["AA", "AA"]

def test_order_by_similarity_three_sequences():
    """
    Testing correct order of sequences based on similarity
    """
    kmer_dict = {
        "A": ["ABC", "BCD"],
        "B": ["ABC", "XYZ"],
        "C": ["FFF"]
    }
    order = order_by_kmer_similarity(kmer_dict)
    assert order[0] in ["A", "B"]  # A and B are the closest pair

def test_build_binary_tree_three_sequences():
    """
    Testing correct creation of binary tree
    """
    terminals, internals = build_binary_tree(["A", "B", "C"])
    assert terminals == {1: "A", 2: "B", 3: "C"}
    assert len(internals) == 2