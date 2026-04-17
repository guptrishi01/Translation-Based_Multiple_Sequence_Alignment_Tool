# Translation-Based Multiple Sequence Alignment Tool

> A Python command-line pipeline that aligns nucleotide sequences at the **protein level** — detecting ORFs, translating to amino acids, pairwise aligning with Needleman-Wunsch, and back-translating to codon-aware DNA alignments. Built from scratch as a final project for Programming II at UNC Charlotte.

[![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white)]()
[![Biopython](https://img.shields.io/badge/Biopython-1.80%2B-006400)]()
[![pytest](https://img.shields.io/badge/tested%20with-pytest-0A9EDC)]()
[![License](https://img.shields.io/badge/License-GPL--3.0-blue)]()

**Author:** Rishi Gupta · [guptrishi01@gmail.com](mailto:guptrishi01@gmail.com)

---

## Table of Contents
- [Why Translation-Based Alignment?](#why-translation-based-alignment)
- [Pipeline](#pipeline)
- [Installation](#installation)
- [Usage](#usage)
- [Input / Output](#input--output)
- [Project Structure](#project-structure)
- [Testing](#testing)
- [AI Usage Statement](#ai-usage-statement)

---

## Why Translation-Based Alignment?

Direct nucleotide alignment ignores the fact that coding sequences evolve under selective pressure at the protein level. Aligning amino acids first respects that biology:

- **Synonymous substitutions don't disrupt alignment.** Two codons encoding the same amino acid collapse to a match in protein space.
- **Frame-shift indels are impossible by construction.** Any gap introduced at the protein level becomes a clean 3-nucleotide gap in DNA space.
- **Conserved regions become visible.** Protein-level conservation is often the signal of biological importance.

This project reimplements the classical Needleman-Wunsch algorithm from first principles while respecting that biological reality.

---

## Pipeline

```
┌──────────────────┐
│  1. Input FASTA  │  Nucleotide sequences
└────────┬─────────┘
         ▼
┌──────────────────┐
│  2. ORF Detection│  Scan all 6 reading frames
└────────┬─────────┘
         ▼
┌──────────────────┐
│  3. Translation  │  Nucleotides → amino acids
└────────┬─────────┘
         ▼
┌──────────────────┐
│  4. K-mer Sort   │  Order sequences for progressive alignment
└────────┬─────────┘
         ▼
┌──────────────────┐
│  5. Pairwise N-W │  Needleman-Wunsch on amino acids
└────────┬─────────┘
         ▼
┌──────────────────┐
│ 6. Back-Translate│  Aligned AAs → codon-aware DNA
└────────┬─────────┘
         ▼
┌──────────────────┐
│  7. Variability  │  Per-codon mismatch/indel rates
│     Analysis     │  CSV + plots
└──────────────────┘
```

---

## Installation

```bash
git clone https://github.com/guptrishi01/Translation-Based_Multiple_Sequence_Alignment_Tool.git
cd Translation-Based_Multiple_Sequence_Alignment_Tool

conda create -n bioinfo_env pandas numpy biopython matplotlib pytest -y
conda activate bioinfo_env
```

---

## Usage

```bash
python3 src/main.py [OPTIONS]
```

### Options

| Flag | Description |
|------|-------------|
| `-h`, `--help` | Show help message |
| `-i`, `--input` | Path to input FASTA file |
| `-o`, `--output` | Path to output FASTA alignment |
| `-r`, `--reference` | Sequence ID to use as reference for mutation-rate analysis |
| `-ma`, `--match` | Score for a match |
| `-mi`, `--mismatch` | Penalty for a mismatch |
| `-in`, `--indel` | Penalty for an insertion/deletion |
| `-t`, `--table` | Translation table (e.g., 1 for standard, 2 for vertebrate mitochondrial) |
| `-k`, `--kmer` | K-mer size for similarity ordering |

### Example

```bash
python3 src/main.py \
    -i data/fly_dna1.fasta \
    -o results/aligned.fasta \
    -ma 1 -mi -1 -in -2 \
    -t 1 -k 3
```

---

## Input / Output

### Input FASTA
- Must contain nucleotide sequences
- ORFs must exist in at least one reading frame

```
>seq1
ATGACCTTGAATG...
>seq2
ATGGGCTTTAG...
```

### Output files

| File | Contents |
|------|----------|
| `results/aligned_amino_acids.fasta` | Amino acid alignment |
| `results/dna_codon_alignment.fasta` | Codon-aware DNA alignment |
| `results/codon_positions.csv` | Per-codon-position mismatch/mutation rate |
| `results/alignment_stats.csv` | Per-position mismatch and indel rates |
| `results/alignment_mutation_rates.png` | Variability visualization plot |

### `codon_positions.csv` columns
`codon_position`, `mismatches`, `total`, `mutation_rate`

### `alignment_stats.csv` columns
`Position`, `MismatchRate`, `IndelRate`

---

## Project Structure

```
Translation-Based_Multiple_Sequence_Alignment_Tool/
├── README.md
├── LICENSE
├── src/
│   └── msaligner/
│       ├── __init__.py
│       ├── orf.py                   # ORF detection in 6 frames
│       ├── kmer_ordering.py         # K-mer similarity + ordering
│       ├── needleman_wunsch.py      # Pairwise alignment
│       ├── back_translate.py        # Codon-aware back-translation
│       └── main.py                  # CLI entry point
├── tests/
│   ├── __init__.py
│   ├── conftest.py
│   ├── orf_test.py
│   ├── kmer_test.py
│   ├── align_test.py
│   ├── test_file.fasta
│   └── input_file.fasta
├── data/
│   └── fly_dna1.fasta
├── docs/
└── results/                         # Created on first run
```

---

## Testing

```bash
pytest tests/ -v
```

Unit tests cover:
- ORF detection (known-ORF inputs, sequences with no ORF, frame edge cases)
- K-mer ordering (symmetry, correct pair identification)
- Alignment correctness (known-alignment fixtures)

---

## AI Usage Statement

I used AI (ChatGPT) **as an assistant, not a code generator.** AI helped with:

- Debugging
- Clarifying biological concepts
- Suggesting efficiency improvements
- Reviewing functions for ordering sequences and back-translation

**No code was copied verbatim.** All final implementations were manually written, tested, and adapted for the project's requirements.

---

## License

GPL-3.0 — see [LICENSE](LICENSE).
