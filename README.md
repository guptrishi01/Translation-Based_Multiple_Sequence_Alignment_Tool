@ -1,153 +1,151 @@
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
This project is a **Translation-Based Multiple Sequence Aligner (MSA)** pipeline written in Python. 
It aligns nucleotides by translating them into amino acids, pairwise aligning at the protein level, 
and then back-translating to codon-aware nucleotide alignments.

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
- Assisting in writing functions for ordering sequences and back translation

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
---

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

---


## Input & Output Formats
### Input FASTA format
a) Must contain nucleotide sequences  
b) ORFs must exist in at least one frame  

Example:
```fasta
>seq1
ATGACCTTGAATG...
>seq2
ATGGGCTTTAG...
```

### Output Files
#### 1. Amino Acid Alignment
`results/aligned_amino_acids.fasta`

#### 2. DNA Codon Alignment
`results/dna_codon_alignment.fasta`

#### 3. Codon Statistics
`results/codon_positions.csv` – Columns:
- codon_position
- mismatches
- total
- mutation_rate

`results/alignment_stats.csv` – Columns:
- Position
- MismatchRate
- IndelRate


#### 4. Codon Variability Plot
`results/alignment_mutation_rates.png` 
- Visualizes mismatch and indel rates across alignment

---

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
|       └── main.py
│  
├── tests/  
│   ├── __init__.py  
│   ├── conftest.py  
│   ├── align_test.py  
│   ├── orf_test.py  
│   ├── kmer_test.py
│   └── test_file.fasta   
│   └── input_file.fasta   
│  
├── data/  
│   └── fly_dna1.fasta 
│  
├── results/  
├── results/  (Created after running main.py)
│   ├── aligned_proteins.fasta  
│   ├── aligned_codons.fasta  
│   ├── codon_stats.csv
