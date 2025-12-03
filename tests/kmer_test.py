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
    assert get_kmers("ABCDE", 3) == ["ABC", "BCD", "CDE"]

def test_get_kmers_short_sequence():
    assert get_kmers("AB", 3) == []

def test_get_kmers_invalid_k():
    with pytest.raises(ValueError):
        get_kmers("ABC", 0)


def test_jaccard_basic():
    k1 = ["AAA", "BBB", "CCC"]
    k2 = ["BBB", "CCC", "DDD"]
    assert jaccard_calculation(k1, k2) == 2/4

def test_jaccard_no_overlap():
    assert jaccard_calculation(["AAA"], ["BBB"]) == 0.0

def test_dataframe_to_dict_basic():
    df = pd.DataFrame({
        "id": ["S1", "S2"],
        "Amino_Acids": ["ABCDE", "AAA"]
    })
    d = dataframe_to_dict(df, k=2)
    assert d["S1"] == ["AB", "BC", "CD", "DE"]
    assert d["S2"] == ["AA", "AA"]

def test_order_by_similarity_three_sequences():
    kmer_dict = {
        "A": ["ABC", "BCD"],
        "B": ["ABC", "XYZ"],
        "C": ["FFF"]
    }
    order = order_by_kmer_similarity(kmer_dict)
    assert order[0] in ["A", "B"]  # A and B are the closest pair

def test_build_binary_tree_three_sequences():
    terminals, internals = build_binary_tree(["A", "B", "C"])
    assert terminals == {1: "A", 2: "B", 3: "C"}
    assert len(internals) == 2