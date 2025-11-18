#!/usr/bin/env python3

"""
This program generates a list of open reading frames obtained from the DNA sequences 
of  a FASTA file.

Author: Rishi Gupta
Contact: rgupta25@charlotte.edu

"""

from Bio import SeqIO
import pandas as pd
import sys

def find_orfs(fasta: str) -> pd.DataFrame:
    """
    Find all open reading frames (ORFs) in a nucleotide sequence on a given strand.

    Args:
        fasta (str): The FASTA file entered by the user on the command line
    Returns:
        pd.DataFrame: A A Pandas DataFrame that contains each ORF, where the ORF is represented as:
            - 'start' (int): Start position
            - 'end' (int): End position
            - 'frame' (int): Reading frame (1, 2, 3)
            - 'sequence' (str): The open reading frame nucleotide sequence.
    """

    """ FASTA dictonary """
    fasta_dict = {
        "start" : [],
        "end" : [],
        "frame" : [],
        "seq" : [],
    }
    
    stop_codon = ["TAG", "TAA", "TGA"]

    """ Iterate through all 3 reading frames """
    try:
        with open(fasta, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                for frame in range(3):
                    """ Go through whole sequence """
                    for i in range(frame, len(record.seq) - 2, 3):
                        if record.seq[i:i + 3] == "ATG":
                            for j in range(i + 3, len(record.seq) - 2, 3):
                                if record.seq[j:j + 3] in stop_codon:
                                    fasta_dict["start"].append(i)
                                    fasta_dict["end"].append(j + 3)
                                    fasta_dict["seq"].append(record.seq[i: j + 3])
                                    fasta_dict["frame"].append(frame)
                                    break # Found ORF, break out of for loop
                        else:
                            continue # No start codon found - continue forward
    except FileNotFoundError:
        sys.stderr.write("FASTA file does not exist")
        raise FileNotFoundError
    
    return pd.DataFrame(fasta_dict) # Returns pandas DataFrame of dictionary