#!/usr/bin/env python3

"""
This program generates a list of open reading frames obtained from the DNA sequences 
of  a FASTA file.

Author: Rishi Gupta
Contact: rgupta25@charlotte.edu

"""

def back_translate_alignment(amino_alignment_dict, df):
    """
    Convert aligned AA sequences back to codon-preserving nucleotide alignments.
    """

    codon_alignment = {}
    
    # Convert all ORFs into lists of STR codons
    nuc_map = {}
    for _, row in df.iterrows():
        seq_id = row["id"]
        orf = str(row["orf"]) 
        codons = [str(orf[i:i+3]) for i in range(0, len(orf), 3)]
        nuc_map[seq_id] = codons

    # Build back-translated alignment
    for seq_id, aa_alignment in amino_alignment_dict.items():

        codons = nuc_map[seq_id]
        codon_idx = 0
        codon_aligned = []

        for aa in aa_alignment:

            if aa == "-":
                codon_aligned.append("---")
            else:
                codon_aligned.append(codons[codon_idx])
                codon_idx += 1

        codon_alignment[seq_id] = "".join(codon_aligned)

    return codon_alignment
