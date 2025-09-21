# Programming II Final Project
# Author: Rishi Gupta

## Project Overview
This project is a **Translation-Based Multiple Sequence Aligner (MSA)** pipeline written in Python.
It aligns nucleotides by translating them into amino acids, pairwise aligning at the protein level, and then back-translating to codon-aware nucleotide alignments.

The pipeline includes these functionalities:
1. **ORF Detection** - Identify open reading frames (ORFs) in nucleotide sequences.  
2. **Translation** – Convert nucleotide sequences into amino acids.  
3. **k-mer Similarity & Sorting** – Sort amino acid sequences for progressive alignment 
4. **Pairwise Alignment** – Utilize the Needleman-Wunsch algorithm for pairwise alignment
5. **Back-Translation** – Convert aligned amino acids into codon-aware nucleotide alignments.  
6. **Codon Statistics & Visualization** – Generate CSV summaries and plots of codon position variability.

## Usage

## Repository structure
msaligner_final_project/
├── README.md
├── LICENSE
├── pseudocode.txt
├── run_pipeline.sh
│
├── src/
│   └── msaligner/
│       ├── __init__.py
│       ├── cli.py
│       ├── io_utils.py
│       ├── orf.py
│       ├── translate.py
│       ├── kmers.py
│       ├── align.py
│       ├── backtranslate.py
│       ├── statistics.py
│       └── plots.py
│
├── tests/
│   ├── test_io_utils.py
│   ├── test_orf.py
│   ├── test_translate.py
│   ├── test_align.py
│   └── test_statistics.py
│
├── data/
│   └── example.fasta
│
├── results/
│   ├── aligned_proteins.fasta
│   ├── aligned_codons.fasta
│   ├── codon_stats.csv
│   └── codon_variability.png
│
└── docs/
    └── FAIR_checklist.md
