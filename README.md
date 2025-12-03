# Multiple Sequence Alignment Pipeline  
### **Final Project – Programming II**  
**Author:** Rishi Gupta  
**Email:** guptrishi01@gmail.com  
**License:** GPL-2.0  

---

## Table of Contents
- Project Overview
- AI Usage Statement
- Installation Instructions
- Usage Examples
- Input & Output Formats
- Repository Structure
- Limitations & Assumptions
- License

---

### Project Overview
This project is a **Translation-Based Multiple Sequence Aligner (MSA)** pipeline written in Python. It aligns nucleotides by translating them into amino acids, pairwise aligning at the protein level, and then back-translating to codon-aware nucleotide alignments.

The pipeline includes these functionalities:
1. **ORF Detection** - Identify open reading frames (ORFs) in nucleotide sequences.  
2. **Translation** – Convert nucleotide sequences into amino acids.  
3. **k-mer Similarity & Sorting** – Sort amino acid sequences for progressive alignment 
4. **Pairwise Alignment** – Utilize the Needleman-Wunsch algorithm for pairwise alignment
5. **Back-Translation** – Convert aligned amino acids into codon-aware nucleotide alignments.  
6. **Codon Statistics & Visualization** – Generate CSV summaries and plots of codon position variability.

---

## AI Usage Statement
I used AI (ChatGPT) **as an assistant, not a code generator**. I used ChatGPT in these ways:

- Debugging
- Clarifying biological concepts
- Suggesting methods to increase code efficiency

**No code was copied verbatim.** All final implementations were manually written, tested, and adapted for the project’s requirements.

---

## Installation
Clone the repository:
```bash
git clone https://github.com/guptrishi01/final-project-Rishi_Gupta.git
cd final-project_Rishi_Gupta/src/msaligner
```

Create conda environment
```bash
conda create --name bioinfo_env pandas numpy biopython matplotlib pytest
conda activate bioinfo_env
```

## Usage

```bash
python3 main.py <options>
```

To see the available options:

```bash
python3 main.py -h
```

Options/Arguments:
- -h, --help (Help message)
- -i, --input (Path to input FASTA file)
- -o, --output (Path to output FASTA alignment file)
- -r, --reference (Sequence ID to use as a reference for mutation-rate analysis)
- -ma, -match (Score for a match)
- -mi, -mismatch (Penalty for a mismatch)
- -in, -indel (Penalty for an insertion/deletion)
- -t, -table (Translation table)
- -k, -kmer (K-mer size)


###


### Project Repository Structure
final-project-Rishi_Gupta/  
├── README.md  
├── LICENSE   
│  
├── src/  
│   └── msaligner/  
│       ├── __init__.py  
│       ├── back_translate.py  
│       ├── kmer_ordering.py  
│       ├── orf.py  
│       ├── needleman_wunsch.py  
│       ├── main.py  
|       └── test.fasta
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
    └── FAIR_checklist

### License

This project operates under the GPL-2.0 license
