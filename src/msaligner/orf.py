#!/usr/bin/env python3

"""
This program generates a list of open reading frames obtained from the DNA sequences 
of  a FASTA file.

Author: Rishi Gupta
Contact: rgupta25@charlotte.edu

"""
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import sys


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
    # Forward direction (5' - 3') is normal transcription and translation, Reverse direction (3' - 5') needs to be reverse complemented to forward
    df['Amino_Acids'] = [ Seq(row.orf).reverse_complement().transcribe().translate(table = user_table) if row.reverse == True 
                   else (Seq(row.orf).transcribe().translate(table = user_table)) 
                   for row in df.itertuples(index=False) ]
    
    return df

def find_orfs(fasta: str) -> pd.DataFrame:
    """
    Find the longest open reading frames (ORFs) for sequences in a FASTA file

    Args:
        fasta (str): The FASTA file entered by the user on the command line
    Returns:
        pd.DataFrame: A Pandas DataFrame that contains the longest ORF of each sequence, where each ORF represented as:
            - 'sequence' (str): The open reading frame nucleotide sequence.
            - 'reverse' (bool): If it's a forward or reverse sequence
            - 'start' (int): Start position
            - 'end' (int): End position
            - 'frame' (int): Reading frame (1, 2, 3)
            - 'orf' (str): The longest open reading frame
    """
    """ List of Longest Open Reading Frame for FASTA Sequences"""
    orf_list = []

    try:
        with open(fasta, "r") as fh:
            orf_list = [find_longest_orf(record.seq).head(1) for record in SeqIO.parse(fh, "fasta")]
    except FileNotFoundError:
        sys.stderr.write("FASTA file was not found")
        raise FileNotFoundError
    except IOError as e:
        sys.stderr.write(f"An IO error has occurred {e}")
        raise IOError
    except Exception as e:
        sys.stderr.write(f"An error has occurred: {e}")
        raise Exception
    return pd.concat(orf_list, ignore_index=True)

def find_longest_orf(seq: str) -> pd.Series:
    """
    Find the longest open reading frames (ORF) in a nucleotide sequence on a given strand.

    Args:
        seq (str): The nucleotide sequence
        fasta_dict (dict): The dictionary that contains the longest ORFs from FASTA file
        reverse (bool): If the nucleotide sequence is forward or reverse
    
    Returns:
        pd.DataFrame : Pandas DataFrame containing ORFs where the longest ORF is at top
    """
    
    stop_codon = ["TAG", "TAA", "TGA"]

    """ FASTA dictonary """
    fasta_dict = {
        "start" : [],
        "end" : [],
        "orf" : [],
        "frame" : [],
        "reverse" : [],
    }

    """ Forward and Reverse Strand """
    for nucleotide in [seq ,seq.reverse_complement()]:
        """ Reading frames 0, 1, 2 """
        for frame in range(3):
            """ Go through whole sequence """
            i = frame
            for i in range(frame, len(nucleotide) - 2, 3):
                """ If a start codon is found - go through rest of sequence for stop codon """
                if nucleotide[i:i + 3] == "ATG":
                    # If another start codon is found in an existing ORF, skip it
                    for j in range(i + 3, len(nucleotide) - 2, 3):
                        """ Add reading frame into FASTA dictionary if stop codon is found """
                        if nucleotide[j:j + 3] in stop_codon:
                            fasta_dict["start"].append(i)
                            fasta_dict["end"].append(j + 3)
                            fasta_dict["orf"].append(nucleotide[i: j + 3])
                            fasta_dict["frame"].append(frame)
                            fasta_dict["reverse"].append(True if nucleotide != seq else False)
                            break # Found ORF, break out of for loop
    pd.DataFrame(fasta_dict).to_csv('output.csv')
    """ Sort Dictionary to where longest ORF is at top """
    return pd.DataFrame(fasta_dict).sort_values(by='orf', key=lambda x: x.str.len(), ascending=False)