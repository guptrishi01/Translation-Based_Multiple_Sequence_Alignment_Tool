#!/usr/bin/env python3

"""
This module provides functions for finding and translating the
longest open reading frames for each sequence of a FASTA file.

It includes functions for parsing a FASTA file, identifying and translating open reading frames,
and returning the longest ORF for each sequence. The goal is to use find these ORFs for each
sequence and pairwise align them via the Needleman Wunsch algorithm.

Functions:
    - find_orfs(fasta, user_table): Finds the longest ORFs 
    - find_longest_orf(seq, id, user_table): Applies specified transformations.
"""

from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import sys

class ORFError(Exception):
    """
    Custom error that is raised due to invalid or missing data
    """
    pass

def find_orfs(fasta: str, user_table : int=1) -> pd.DataFrame:
    """
    Find the longest open reading frames (ORFs) for sequences in a FASTA file

    Args:
        fasta (str): The FASTA file entered by the user on the command line
        user_table (int): The genetic code to use for translation (default is 1)
    Returns:
        pd.DataFrame: A Pandas DataFrame that contains the longest ORF of each sequence, where each ORF represented as:
            - 'id' (str) : Sequence ID
            - 'start' (int): Start position
            - 'end' (int): End position
            - 'sequence' (str): The DNA sequence that was scanned
            - 'orf' (str): The longest open reading frame
            - 'frame' (int): Reading frame (1, 2, 3)
            - 'Amino_Acids' (Seq): The translated longest open reading frame
    """
    orf_list = []

    # Reading in FASTA file + Error Handling

    try:
        with open(fasta, "r") as fh:
            orf_list = [find_longest_orf(record.seq, record.id, user_table).head(1) for record in SeqIO.parse(fh, "fasta")]
    
    except FileNotFoundError:
        sys.stderr.write("FASTA file was not found")
        raise FileNotFoundError("FASTA file was not found")
    
    except IOError as e:
        sys.stderr.write(f"An IO error has occurred {e}")
        raise IOError(f"An IO error has occurred {e}")
    
    except ORFError:
        sys.stderr.write(f"Failed to parse FASTA file")
        raise Exception(f"Failed to parse FASTA file")
    
    if not orf_list:
        sys.stderr.write("FASTA file contained no sequences")
        raise ORFError("FASTA file contained no sequences")
    
    return pd.concat(orf_list, ignore_index=True)

def find_longest_orf(seq: str, id : str, user_table : int=1) -> pd.DataFrame:
    """
    Find the longest open reading frames (ORF) in a nucleotide sequence on a given strand.

    Args:
        seq (str): The nucleotide sequence
        id (str): The nucleotide sequence ID
        user_table (int): The genetic code to use for translation (default is 1)
    
    Returns:
        pd.DataFrame : Pandas DataFrame containing ORFs where the longest ORF is at top
    """
    bases = ['A', 'T', 'C', 'G']

    for read in seq:
        if read not in bases:
            sys.stderr.write("Invalid base in DNA sequence")
            raise ORFError("Invalid base in DNA sequence")

    if seq == "":
        sys.stderr.write("Sequence from FASTA file is empty")
        raise ORFError("Sequence from FASTA file is empty")
    
    stop_codon = ["TAG", "TAA", "TGA"]

    # FASTA dictionary
    fasta_dict = {
        "id" : [],
        "start" : [],
        "end" : [],
        "seq" : [],
        "orf" : [],
        "frame" : [],
        "Amino_Acids" : []
    }

    # Go through all 6 reading frames
    for nucleotide in [seq ,seq.reverse_complement()]:
        for frame in range(3):
            i = frame
            for i in range(frame, len(nucleotide) - 2, 3):
                # Go through rest of sequence once start codon is found
                if nucleotide[i:i + 3] == "ATG":
                    # If another start codon is found in an existing ORF, skip it
                    for j in range(i + 3, len(nucleotide) - 2, 3):
                        # Add ORF if stop codon is found
                        if nucleotide[j:j + 3] in stop_codon:
                            fasta_dict["id"].append(id)
                            fasta_dict["start"].append(i)
                            fasta_dict["end"].append(j + 3)
                            fasta_dict["seq"].append(nucleotide)
                            fasta_dict["orf"].append(nucleotide[i: j + 3])
                            fasta_dict["frame"].append(frame)
                            fasta_dict["Amino_Acids"].append(Seq(nucleotide[i: j + 3]).translate(table=user_table))
                            break
    if not fasta_dict:
        sys.stderr.write("No ORF was found in sequence")
        raise ORFError("No ORF was found in sequence")
    # Sort dictionary where longest ORF is at top
    return pd.DataFrame(fasta_dict).sort_values(by='orf', key=lambda x: x.str.len(), ascending=False)