#!/usr/bin/env python3

import unittest
from Bio import SeqIO
import pandas as pd
import sys

def find_orfs(fasta: str) -> pd.DataFrame:
    orf_list = []

    try:
        with open(fasta, "r") as fh:
            for record in SeqIO.parse(fh, "fasta"):
                orf_list.append(find_longest_orf(record.seq))
    except FileNotFoundError:
        sys.stderr.write("FASTA file was not found")
        raise FileNotFoundError
    except IOError as e:
        sys.stderr.write(f"An IO error has occurred {e}")
        raise IOError
    except Exception as e:
        sys.stderr.write(f"An error has occurred: {e}")
        raise Exception
    print(len(orf_list))
    return pd.concat(orf_list, axis=1)

def find_longest_orf(seq: str) -> pd.Series:
    stop_codon = ["TAG", "TAA", "TGA"]
    fasta_dict = {
        "seq" : [],
        "start" : [],
        "end" : [],
        "orf" : [],
        "frame" : [],
        "reverse" : [],
    }

    for nucleotide in [seq, seq[::-1]]:
        for frame in range(3):
            for i in range(frame, len(nucleotide) - 2, 3):
                if nucleotide[i:i + 3] == "ATG":
                    for j in range(i + 3, len(nucleotide) - 2, 3):
                        if nucleotide[j:j + 3] in stop_codon:
                            fasta_dict["seq"].append(nucleotide)
                            fasta_dict["start"].append(i)
                            fasta_dict["end"].append(j + 3)
                            fasta_dict["orf"].append(nucleotide[i: j + 3])
                            fasta_dict["frame"].append(frame)
                            fasta_dict["reverse"].append(True if nucleotide != seq else False)
                            break
                        else:
                            continue
    
    """ Sort Dictionary to grab longest ORF """
    return pd.DataFrame(fasta_dict).sort_values(by='orf', key=lambda x: x.str.len(), ascending=False).iloc(0)

class FindORFTestCase(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()