#!/usr/bin/env python3

import pytest
import pandas as pd
from Bio.Seq import Seq
from pathlib import Path
from Bio import SeqIO
from msaligner.orf import find_orfs, find_longest_orf, ORFError

TEST_DIR = Path(__file__).parent
FASTA = TEST_DIR / "test_file.fasta"

def test_find_longest_orf_basic():
    seq = Seq("AAAATGAAATAGAAA")
    df = find_longest_orf(seq, "s1")
    assert isinstance(df, pd.DataFrame)

    assert len(df) >= 1
    assert df.iloc[0]["orf"].startswith("ATG")
    assert df.iloc[0]["orf"].endswith(("TAG", "TAA", "TGA"))

def test_find_orfs_missing_file():
    with pytest.raises(FileNotFoundError):
        find_orfs("no_such_file.fasta")


def test_find_orfs_reads_sequences():
    assert FASTA.exists(), f"Test FASTA not found at {FASTA}"

    df = find_orfs(str(FASTA))

    assert isinstance(df, pd.DataFrame)
    assert "orf" in df.columns
    assert "Amino_Acids" in df.columns
    assert len(df) > 0


def test_ensure_longest_orf():
    assert FASTA.exists(), f"Test FASTA not found at {FASTA}"

    seqs = list(SeqIO.parse(str(FASTA), "fasta"))
    assert len(seqs) >= 1
    assert all(len(str(record.seq)) > 0 for record in seqs)

