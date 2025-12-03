#!/usr/bin/env python3

import pytest
import pandas as pd
from Bio.Seq import Seq
from msaligner.orf import find_orfs, find_longest_orf, ORFError


def test_find_longest_orf_basic():
    seq = "AAAATGAAATAGAAA"
    df = find_longest_orf(seq, "s1")
    assert not df.empty
    assert df.iloc[0]["orf"] == "ATGAAATAG"


def test_find_longest_orf_no_orf():
    with pytest.raises(ORFError):
        seq = "AAAAAAAAAAAAAA"
        df = find_longest_orf(seq, "s1")
        assert df.empty


def test_find_orfs_missing_file():
    with pytest.raises(FileNotFoundError):
        find_orfs("no_such_file.fasta")


def test_find_orfs_reads_sequences():
    df = find_orfs("test.fasta")
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 3


def test_ensure_longest_orf():
    df_all = []

    from Bio import SeqIO
    for record in SeqIO.parse("test.fasta", "fasta"):
        df = find_longest_orf(record.seq, record.id)
        if not df.empty:
            df_all.append((record.id, df))

    df_final = find_orfs("test.fasta")

    for seq_id, df in df_all:
        longest_true = df.iloc[0]["orf"]
        longest_tested = df_final[df_final["id"] == seq_id].iloc[0]["orf"]
        assert len(longest_tested) == len(longest_true)

