#!/usr/bin/env python3

"""
This program imports the other files in this
subdirectory and utilizes the command line interface
to take in an input FASTA file and perform the appropriate
calculates and visualizations.

Author: Rishi Gupta
Contact: rgupta25@charlotte.edu
"""

import pandas as pd

import argparse, sys
from orf import find_orfs
from translate_to_AA import protein_translate


# def get_args():
# 	"""Parse command-line arguments for grid dimensions."""
# 	parser = argparse.ArgumentParser(
# 		description="Draw an ASCII grid with customizable dimensions."
# 	)
	
# 	parser.add_argument(
# 		"--rows",
# 		type=int,
# 		default=2,
# 		help="Number of rows in the grid (default: 2)"
# 	)

# 	return parser.parse_args()

def main():
    """
    Utilizes functions from other programs detect ORFs, translate sequences
    to amino acids, perform pairwise alignment, and back-translate alignments and
    generate codon position statistics and visualizations.
    """
    test_fasta = "fly_dna1.fasta"
    orfs = find_orfs(test_fasta)
    print(len(orfs))
    print(protein_translate(orfs))
    
if __name__ == "__main__":
    main()