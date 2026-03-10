#!/usr/bin/env python3


"""
This module provides an end-to-end workflow for performing multiple
sequence alignment from raw nucleotide FASTA input to aligned DNA
open reading frames.


Workflow is encapsulated in the MSAligner class, which manages:
	- ORF Detection
	- Amino-Acid Translation
	- K-mer detection and Guide Tree Construction
	- Pairwise Needleman-Wunsch Alignment
	- Consensus and Gap Patching
	- Codon Back-Translation
	- Output Writing and Visualizations

Author: Rishi Gupta
Contact: guptrishi01@gmail.com
"""

import pandas as pd
import numpy as np
import argparse, sys
import os
from typing import List, Dict
from pathlib import Path
import matplotlib.pyplot as plt

from orf import find_orfs
from kmer_ordering import dataframe_to_dict
from kmer_ordering import order_by_kmer_similarity
from kmer_ordering import build_binary_tree
from needleman_wunsch import init_mat
from needleman_wunsch import fill_matrix
from needleman_wunsch import trace_matrix

# Project Root and Output Directory Path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR.mkdir(exist_ok=True)


def parse_args() -> argparse.Namespace:
	"""Parse command-line arguments for the program."""
	parser = argparse.ArgumentParser(
		description="Multiple sequence alignment tool."
	)
	parser.add_argument(
		"-i", "--input",
		default="input_file.fasta",
		help="Path to input FASTA file."
	)
	parser.add_argument(
		"-o", "--output",
		default="aligned_amino_acids.fasta",
		help="Path to output FASTA alignment file."
	)
	parser.add_argument(
		"-r", "--reference",
		default="Seq1",
		help="Sequence ID to use as reference for mutation-rate analysis."
	)
	parser.add_argument("-ma", "--match", default=1, type=int, help="Score for a match.")
	parser.add_argument("-mi", "--mismatch", default=-1, type=int, help="Penalty for a mismatch.")
	parser.add_argument("-in", "--indel", default=-1, type=int, help="Penalty for an insertion/deletion.")
	parser.add_argument("-t", "--table", default=1, type=int, help="Translation table.")
	parser.add_argument("-k", "--kmer", default=3, type=int, help="K-mer size.")
	return parser.parse_args()



class MSAligner:
	def __init__(self, fasta_file: str, aln_file: str, 
				 ref : str, match : int, mismatch : int,
				 indel : int, ttable : int, ksize : int) -> None:
		"""
        Initialize settings and parameters for an MSAligner instance.

        Args: 
			fasta_file (str) : Path to input nucleotide FASTA file.
			aln_file (str) : Path where the final aligned FASTA will be written.
			ref (str) : Sequence ID to treat as reference in downstream mutation analysis.
			match (int) : Alignment score for amino-acid matches.
			mismatch (int) : Alignment penalty for amino-acid mismatches.
			indel (int) : Gap insertion/deletion penalty.
			ttable (int) : NCBI translation table to use for ORF translation.
			ksize (int) : K-mer size used for similarity ordering and guide-tree construction.
        """
		self.fasta = Path(fasta_file)
		self.out = Path(aln_file)
		self.reference = ref
		self.match = match
		self.mismatch = mismatch
		self.indel = indel
		self.table = ttable
		self.kmer = ksize
		
	
	def detect_orfs(self):
		"""
		Detect open reading frames in input FASTA

		"""

		if not self.fasta.exists():
			raise FileNotFoundError(f"FASTA file not found: {self.fasta}")
		self.df = find_orfs(str(self.fasta), self.table)

	def find_kmer_translate(self):
		"""
		Convert ORF data into k-mer dictionary and build guide tree
		
		"""
		self.dict = dataframe_to_dict(self.df, self.kmer)
		self.order = order_by_kmer_similarity(self.dict)
		self.terminals, self.nodes = build_binary_tree(self.order)
		self.sequences = {
		seq_id: np.array([{aa} for aa in self.dict[seq_id]], dtype=object)
		for seq_id in self.dict
		}

	def matrix_initialization(self, seq1, seq2):
		"""
		Initialize Needleman-Wunsch scoring matrix
		
		Args:
		seq1 : Array of character sets for sequence 1
		seq2 : Array of character sets for sequence 2
		"""

		self.seq1 = seq1
		self.seq2 = seq2
		self.mat = init_mat(seq1, seq2)
	
	def fill(self) -> None:
		"""
		Fill scoring matrix using Needleman-Wunsch algorithm
		
		"""
		self.mat = fill_matrix(
			matrix=self.mat,
			seq1_array=self.seq1,
			seq2_array=self.seq2,
			match=self.match,
			mismatch=self.mismatch,
			indel=self.indel
		)

	def trace(self) -> None:
		"""
		Perform traceback through the filled DP matrix
		to reconstruct the aligned sequences.
		"""
		self.aln = trace_matrix(
			matrix=self.mat,
			seq1_array=self.seq1,
			seq2_array=self.seq2
		)

	def update_sequences(self, id1: str, id2: str) -> None:
		"""
		Update sequences after alignment
		
		id1 : First sequence
		id2: Second sequence
		"""
		self.sequences[id1] = np.array(self.aln[0], dtype=object)
		self.sequences[id2] = np.array(self.aln[1], dtype=object)

	def create_consensus(self, id1: str, id2: str, node: int) -> None:
		"""
		Create consensus sequence for an internal node
		
		id1 : First child sequence
		id2: Second child sequence
		node : Internal node identifier
		"""
		seq1, seq2 = self.aln
		consensus = [seq1[i].union(seq2[i]) for i in range(len(seq1))]
		self.sequences[node] = np.array(consensus, dtype=object)
	
	def patch(self, history: list) -> None:
		"""
		Insert new gaps into previously aligned sequences
		
		history : List of sequence/node IDs that have already been aligned
		"""
		gap_positions = sorted([i for i, char in enumerate(self.aln[0]) if char == {'-'}], reverse=True)
		for seq_id in history:
			seq = self.sequences[seq_id]
			seq_len = len(seq)
			trailing = []
			# Insert gaps at the same positions
			for pos in gap_positions:
				if pos > seq_len:
					trailing.append(pos)
				else:
					seq = np.insert(seq, pos, {'-'})

			# If gaps would occur past the end, append them
			if trailing:
				seq = np.append(seq, [{' - '} for _ in trailing])

			self.sequences[seq_id] = seq
			
	def back_translate_alignment(self):
		"""
		Convert amino acid alignment into codon-level alignment
		
		"""

		# Dictionary that will contain codons of DNA ORF for each sequence
		nuc_map = {}

		# Iterate through DataFrame containing DNA sequence ID and ORF 
		for _, row in self.df.iterrows():
			seq_id = row["id"]
			orf = str(row["orf"])                      # ensure it's str
			codons = [orf[i:i+3] for i in range(0, len(orf), 3)]
			nuc_map[seq_id] = codons

		
		codon_alignment = {}

		# Go through each amino acid alignment
		for seq_id, aa_aln in self.sequences.items():

			# Only terminals have codons
			if seq_id not in nuc_map:
				continue

			codons = nuc_map[seq_id]
			codon_idx = 0
			codon_aligned = []

			for aa_set in aa_aln:

				# Convert set to a usable symbol:
				aa = next(iter(aa_set))

				# If gap, produce a 3 nucleotide gap codon (---)
				if aa == '-':
					codon_aligned.append('---')
				else:
					# If codons are still available
					if codon_idx < len(codons):
						codon_aligned.append(codons[codon_idx])
						codon_idx += 1
					else:
						# Fill  with gap codons if alignment is longer than ORF
						codon_aligned.append("---")

			# Dictionary mapping sequence ID to codon-aligned nucleotide
			codon_alignment[seq_id] = "".join(codon_aligned)

		# Store output
		self.codon_alignment = codon_alignment
	
	def codon_position_statistics(self):
		"""
		Calculate mutation rates by codon position (1, 2, 3) using
		codon-aligned nucleotide sequences
		"""

		ref_id = self.reference

		# Error handling - if user input reference is not found in alignment
		if ref_id not in self.codon_alignment:
			raise ValueError(f"Reference {self.reference} not in codon alignment")

		ref = self.codon_alignment[ref_id]
		seqs = {seq_id : seq for seq_id, seq in self.codon_alignment.items() if seq_id != ref_id}
		pos_counts = {1: {'mismatch': 0, 'total': 0},
                  2: {'mismatch': 0, 'total': 0},
                  3: {'mismatch': 0, 'total': 0}}

		# walk codons
		for i in range(0, len(ref), 3):
			r1, r2, r3 = ref[i:i+3]

			for seq_id, seq in seqs.items():
				c1, c2, c3 = seq[i:i+3]

				# pos1
				if r1 != '-' and c1 != '-':
					pos_counts[1]['total'] += 1
					if r1 != c1:
						pos_counts[1]['mismatch'] += 1

				# pos2
				if r2 != '-' and c2 != '-':
					pos_counts[2]['total'] += 1
					if r2 != c2:
						pos_counts[2]['mismatch'] += 1

				# pos3
				if r3 != '-' and c3 != '-':
					pos_counts[3]['total'] += 1
					if r3 != c3:
						pos_counts[3]['mismatch'] += 1

		# Build DataFrame
		df = pd.DataFrame({
			"codon_position": [1, 2, 3],
			"mismatches": [pos_counts[i]['mismatch'] for i in (1, 2, 3)],
			"total": [pos_counts[i]['total'] for i in (1, 2, 3)]
		})

		df["mutation_rate"] = df["mismatches"] / df["total"]

		self.codon_stats_df = df
		csv_path = OUTPUT_DIR / "codon_position_stats.csv"
		df.to_csv(csv_path, index=False)
		sys.stdout.write("Saved codon_position_stats.csv")

	def alignment_statistics(self):
		"""
		Compute mismatch and indel rates at every alignment position relative to 
		the reference sequence. Works on the codon-aligned nucleotide sequences.
		"""

		data = self.codon_alignment
		ref_id = self.reference

		ref_seq = data[ref_id]
		seq_ids = [sid for sid in data if sid != ref_id]

		n_positions = len(ref_seq)
		n_seqs = len(seq_ids)


		# Counters
		mismatch_counts = [0] * n_positions
		indel_counts = [0] * n_positions

		# Compare each sequence to reference
		for sid in seq_ids:
			alt = data[sid]

			for i, (r, b) in enumerate(zip(ref_seq, alt)):

				# both gaps → ignore
				if r == "-" and b == "-":
					continue

				# indel event
				if r == "-" or b == "-":
					indel_counts[i] += 1

				# mismatch event
				elif r != b:
					mismatch_counts[i] += 1

		# Per-position rates
		mismatch_rates = [m / n_seqs for m in mismatch_counts]
		indel_rates    = [g / n_seqs for g in indel_counts]

		# DataFrame output
		df = pd.DataFrame({
			"Position": range(1, n_positions + 1),
			"MismatchRate": mismatch_rates,
			"IndelRate": indel_rates
		})

		self.aln_stats_df = df

		# Compute alignment-wide averages
		self.avg_mismatch = sum(mismatch_counts) / (n_positions * n_seqs)
		self.avg_indel    = sum(indel_counts) / (n_positions * n_seqs)

		# Save results
		out_csv = OUTPUT_DIR / "alignment_stats.csv"
		df.to_csv(out_csv, index=False)

		return df, self.avg_mismatch, self.avg_indel
	
	def plot_alignment_stats(self):
		df = self.aln_stats_df

		plt.figure(figsize=(12,4))
		plt.plot(df["Position"], df["MismatchRate"], label="Mismatch Rate", linewidth=1.2)
		plt.plot(df["Position"], df["IndelRate"], label="Indel Rate", linewidth=1.2)

		plt.xlabel("Alignment Position")
		plt.ylabel("Rate")
		plt.title("Mismatch and Indel Rates Across Alignment")
		plt.legend()
		plt.tight_layout()

		out_png = OUTPUT_DIR / "alignment_mutation_rates.png"
		plt.savefig(out_png, dpi=140)
		plt.show()
		plt.close()

	
	def _tree_root(self):
		"""Return the root node identifier for the guide tree."""
		if not self.nodes:
			# No internal nodes: only one terminal
			return list(self.terminals.values())[0]
		return max(self.nodes.keys())

	def format_guide_tree_newick(self) -> str:
		"""
		Return a Newick-like string representation of the guide tree.
		Leaves are sequence IDs. Internal nodes are labeled N<id>.
		"""
		def rec(node) -> str:
			# Leaves are sequence labels (strings)
			if isinstance(node, str):
				return node
			# Internal nodes are ints pointing into self.nodes
			if isinstance(node, int) and node in self.nodes:
				left, right = self.nodes[node]
				return f"({rec(left)},{rec(right)})N{node}"
			# Fallback (shouldn't happen if the tree is consistent)
			return str(node)

		root = self._tree_root()
		return rec(root) + ";"

	def print_guide_tree(self) -> None:
		"""
		Print ordering + terminals + internal nodes + Newick-like tree to stdout.
		Call this BEFORE alignment output.
		"""
		sys.stdout.write("Guide Tree (k-mer similarity progressive tree)\n")
		sys.stdout.write(f"Ordering: {self.order}\n\n")

		sys.stdout.write("Terminal nodes (id -> sequence):\n")
		for tid in sorted(self.terminals.keys()):
			sys.stdout.write(f"  {tid}: {self.terminals[tid]}\n")

		sys.stdout.write("\nInternal nodes (node_id -> [left, right]):\n")
		for nid in sorted(self.nodes.keys()):
			sys.stdout.write(f"  N{nid}: {self.nodes[nid]}\n")

		sys.stdout.write("\nNewick-like representation:\n")
		sys.stdout.write(self.format_guide_tree_newick() + "\n\n")


	def write_fa(self) -> None:
		"""
		Write final aligned amino-acid and DNA codon-aligned sequences to FASTA format,
		and print them to stdout in two separate sections (AA first, then codon).
		"""
		dna_path = OUTPUT_DIR / "dna_codon_alignment.fasta"
		aa_path = self.out

		# Grab Sequence IDs and Alignments (terminals only)
		seqids = [self.terminals[key] for key in self.terminals]
		aa_seqs = [''.join(next(iter(s)) for s in self.sequences[seqid]) for seqid in seqids]

		if len(seqids) != len(aa_seqs):
			raise ValueError("IDs and sequences must have the same length.")

		# Clean AA sequences once
		aa_clean = {sid: seq.strip().upper() for sid, seq in zip(seqids, aa_seqs)}

		# --------------------------
		# 1) Write AA alignment file
		# --------------------------
		with open(aa_path, "w", encoding="utf-8") as f:
			for sid in seqids:
				f.write(f">{sid}\n{aa_clean[sid]}\n")

		# --------------------------
		# 2) Print AA alignment section
		# --------------------------
		sys.stdout.write("Amino Acid Alignment\n")
		for sid in seqids:
			sys.stdout.write(f">{sid}\n{aa_clean[sid]}\n")

		sys.stdout.write(f"\nAlignment saved as {aa_path}\n\n\n")

		# -----------------------------------
		# 3) Write Codon alignment file (only codon records)
		# -----------------------------------
		with open(dna_path, "w", encoding="utf-8") as f:
			for sid in seqids:
				codon_seq = self.codon_alignment.get(sid, "").strip().upper()
				f.write(f">{sid}_codon\n{codon_seq}\n")

		# -----------------------------------
		# 4) Print Codon alignment section
		# -----------------------------------
		sys.stdout.write("Codon Alignment\n")
		for sid in seqids:
			codon_seq = self.codon_alignment.get(sid, "").strip().upper()
			sys.stdout.write(f">{sid}_codon\n{codon_seq}\n")

		sys.stdout.write(f"\nAlignment saved as {dna_path}\n//\n")

	

def main():
	"""
	Utilizes functions from other programs detect ORFs, translate sequences
	to amino acids, perform pairwise alignment, and back-translate alignments and
	generate codon position statistics and visualizations.
	"""
	args = parse_args()

	# Resolve input path
	input_path = Path(args.input)

	if not input_path.is_absolute():
		input_path = DATA_DIR / input_path
	
	if not input_path.exists():
		raise FileNotFoundError(f"Input FASTA not found: {input_path}")
	
	# Resolve output path
	output_path = OUTPUT_DIR / args.output

	# Create class
	msaligner = MSAligner(
		fasta_file=input_path,
		aln_file=output_path,
		ref=args.reference,
		match=args.match,
		mismatch=args.mismatch,
		indel=args.indel,
		ttable=args.table,
		ksize=args.kmer
	)

	# Utilize class functions
	history: List[str] = []
	msaligner.detect_orfs()
	msaligner.find_kmer_translate()
	msaligner.print_guide_tree()

	# Sequentially align pairs according to the progressive tree
	for internal_node in msaligner.nodes:
		head1, head2 = msaligner.nodes[internal_node]
		seq1, seq2 = msaligner.sequences[head1], msaligner.sequences[head2]

		msaligner.matrix_initialization(seq1, seq2)
		msaligner.fill()
		msaligner.trace()
		msaligner.update_sequences(head1, head2)
		msaligner.create_consensus(head1, head2, internal_node)
		msaligner.patch(history)

		# Track all sequences that have been processed
		for name in (head1, head2):
			if name not in history:
				history.append(name)

	msaligner.back_translate_alignment()
	# Write final alignment to disk
	msaligner.write_fa()
	msaligner.codon_position_statistics()
	msaligner.alignment_statistics()
	msaligner.plot_alignment_stats()
	
	
if __name__ == "__main__":
	main()