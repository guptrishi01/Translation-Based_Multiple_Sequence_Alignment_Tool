#!/usr/bin/env python3

"""
This program generates a list of open reading frames obtained from the DNA sequences 
of  a FASTA file.

Author: Rishi Gupta
Contact: rgupta25@charlotte.edu

"""

from Bio import SeqIO
from Bio.SeqUtils import orf

def find_orfs(fasta: str, min_length: int = 100) -> list[dict]:
    """
    Find all open reading frames (ORFs) in a nucleotide sequence on a given strand.

    Args:
        fasta (str): The FASTA file entered by the user on the command line
        min_length (int): The minimum number of nucleotides in the ORF (the default is 100).

    Returns:
        list[dict]: A list of ORFs, where each ORF is represented as a dictionary containing:
            - 'start' (int): Start position (0-based index).
            - 'end' (int): End position (exclusive).
            - 'frame' (int): Reading frame (0, 1, or 2).
            - 'strand' (str): '+' or '-'.
            - 'sequence' (str): The ORF nucleotide sequence.
    """
    # Create 
    # Use try-except statement to use SeqIO.parse() to read in FASTA file and create iterator
    try:
        sequence_records = SeqIO.parse(fasta, "fasta")
    except FileNotFoundError:
        raise("Error: The file entered could not be found")
    
    # Go through iterator using a for loop and find orf
    # Create entry and add it into a list containing start, end, sequence


