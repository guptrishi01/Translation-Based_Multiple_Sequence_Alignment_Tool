#!/usr/bin/env python3

"""
This program translates nucleotide sequences into amino acids

Author: Rishi Gupta
Contact: rgupta25@charlotte.edu

"""

from Bio.Seq import Seq
import pandas as pd

def protein_translate(df: pd.DataFrame, user_table=1) -> pd.DataFrame:
    """
    Read in DataFrame containing Open Reading Frames and transcribe
    nucleotide sequences to RNA, then RNA sequences to amino acids

    Args:
        df (pd.DataFrame): The Pandas DataFrame containing Open Reading Frames and their DNA sequences
        user_table (int): The user's choice for translation table (default choice is 1)
    
    Returns:
        pd.DataFrame: An updated Pandas DataFrame with a new column
            - 'Amino_Acids' (str): Sequence of amino acids translated from DNA sequence
    """

    amino_acids = []
    for sequence in df['seq']:
        rna_seq = Seq(sequence).transcribe()
        amino_acids.append(rna_seq.translate(table = user_table))
    
    df['Amino_Acids'] = amino_acids

    return df