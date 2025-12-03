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
I used AI (ChatGPT) **as an assistant, not a code generator**. AI helped in these ways:

- Reviewing code structure and identifying logical inconsistencies.
- Suggesting docstring improvements and clarifying function behavior.
- Helping debug algorithmic issues in the alignment and back-translation steps.

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
pip install -r requirements.txt
```

Make the pipeline executable:
```bash
chmod +x run_pipeline.sh
```

Run the full pipeline:
```bash
./run_pipeline.sh data/example.fasta
```

---


### Project Repository Structure
final-project-Rishi_Gupta/  
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
    └── FAIR_checklist

### License

This project operates under the GPL-2.0 license
