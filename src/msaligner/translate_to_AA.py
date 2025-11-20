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

    """ Transcription and Translation of DNA into Amino Acids based on ordering of sequence """
    df['Amino_Acids'] = [ Seq(row.orf).transcribe().translate(table = user_table) if row.reverse == True 
                   else (Seq(row.orf).reverse_complement().transcribe().translate(table = user_table)) 
                   for row in df.itertuples(index=False) ]
    
    return df