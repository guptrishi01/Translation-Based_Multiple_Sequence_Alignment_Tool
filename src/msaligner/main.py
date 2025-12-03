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
import numpy as np
import argparse, sys
from typing import List, Dict
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap
from orf import find_orfs
from kmer_ordering import dataframe_to_dict
from kmer_ordering import order_by_kmer_similarity
from kmer_ordering import build_binary_tree
from needleman_wunsch import init_mat
from needleman_wunsch import fill_matrix
from needleman_wunsch import trace_matrix

def parse_args() -> argparse.Namespace:
	"""Parse command-line arguments for the program."""
	parser = argparse.ArgumentParser(
		description="Multiple sequence alignment tool."
	)
	parser.add_argument(
		"-i", "--input",
		default="test.fasta",
		help="Path to input FASTA file."
	)
	parser.add_argument(
		"-o", "--output",
		default="alignment.fasta",
		help="Path to output FASTA alignment file."
	)
	parser.add_argument(
		"-r", "--reference",
		default="s1",
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
		self.fasta = fasta_file
		self.out = aln_file
		self.reference = ref
		self.match = match
		self.mismatch = mismatch
		self.indel = indel
		self.table = ttable
		self.kmer = ksize
		
	
	def detect_orfs(self):
		self.df = find_orfs(self.fasta, self.table)
	
	def count_kmers(self):
		overall = {}      # per-seq kmer counts
		global_counts = {}  # across all sequences

		for seq_id, kmers in self.dict.items():
			per_seq = {}

			# Count kmers for this sequence
			for k in kmers:
				per_seq[k] = per_seq.get(k, 0) + 1
				global_counts[k] = global_counts.get(k, 0) + 1

			# Sort descending by frequency
			sorted_seq = sorted(per_seq.items(), key=lambda x: x[1], reverse=True)
			overall[seq_id] = sorted_seq

		# Save results
		self.all_kmers = {
			"per_sequence": overall,
			"global": sorted(global_counts.items(), key=lambda x: x[1], reverse=True)
		}



	def build_kmer_matrix(self, top_k=8):
		per_seq = self.all_kmers["per_sequence"]
		global_sorted = self.all_kmers["global"]

		# Select top K kmers across all sequences
		top_kmers = [k for (k, _) in global_sorted[:top_k]]

		matrix = []
		row_names = []

		for seq_id, kmer_list in per_seq.items():
			# Convert sequence list of (kmer, count) → dict
			d = dict(kmer_list)
			row = [d.get(k, 0) for k in top_kmers]
			matrix.append(row)
			row_names.append(seq_id)
		return row_names, top_kmers, np.array(matrix)

	def plot_kmer_heatmap(self, top_k=8):
		self.row_names, self.col_names, M = self.build_kmer_matrix(top_k)
		plt.figure(figsize=(12, 6))
		
		cmap = LinearSegmentedColormap.from_list("blue_soft-orange",["#C8A2C8", "#7a5195", "#4B0082"])

		norm = plt.Normalize(vmin=0, vmax=M.max())
		plt.imshow(M, aspect='auto', cmap=cmap, norm=norm)


		cbar=plt.colorbar(label='k-mer frequency')

	# dynamic ticks based on actual data range
		ticks = list(range(0, int(M.max()) + 1))
		cbar.set_ticks(ticks)
		cbar.set_ticklabels([str(t) for t in ticks])

		

		plt.xticks(range(len(self.col_names)), self.col_names, rotation=90)
		plt.yticks(range(len(self.row_names)), self.row_names)
		plt.title(f"Top {top_k} k-mer Frequency Heatmap")
		plt.tight_layout()
		plt.show()



	def find_kmer_translate(self):
		self.dict = dataframe_to_dict(self.df, self.kmer)
		self.order = order_by_kmer_similarity(self.dict)
		self.terminals, self.nodes = build_binary_tree(self.order)
		self.sequences = {
		seq_id: np.array([{aa} for aa in self.dict[seq_id]], dtype=object)
		for seq_id in self.dict
		}

	def matrix_initialization(self, seq1, seq2):
		self.seq1 = seq1
		self.seq2 = seq2
		self.mat = init_mat(seq1, seq2)
	
	def fill(self) -> None:
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
		self.sequences[id1] = np.array(self.aln[0], dtype=object)
		self.sequences[id2] = np.array(self.aln[1], dtype=object)

	def create_consensus(self, id1: str, id2: str, node: int) -> None:
		seq1, seq2 = self.aln
		consensus = [seq1[i].union(seq2[i]) for i in range(len(seq1))]
		self.sequences[node] = np.array(consensus, dtype=object)
	
	def patch(self, history: list) -> None:
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

		# --------------------------
		# 1. Build codon lookup map
		# --------------------------
		nuc_map = {}

		for _, row in self.df.iterrows():
			seq_id = row["id"]
			orf = str(row["orf"])                      # ensure it's str
			codons = [orf[i:i+3] for i in range(0, len(orf), 3)]
			nuc_map[seq_id] = codons

		# --------------------------
		# 2. Walk AA alignment
		# --------------------------
		codon_alignment = {}

		for seq_id, aa_aln in self.sequences.items():

			# Only terminals have codons
			if seq_id not in nuc_map:
				continue

			codons = nuc_map[seq_id]
			codon_idx = 0
			codon_aligned = []

			for aa_set in aa_aln:

				# Convert set to a usable symbol:
				# e.g. {'A'} → 'A', {'-'} → '-'
				aa = next(iter(aa_set))

				if aa == '-':
					# Amino-acid gap → 3-nt codon gap
					codon_aligned.append('---')
				else:
					# Use next real codon
					codon_aligned.append(codons[codon_idx])
					codon_idx += 1

			codon_alignment[seq_id] = "".join(codon_aligned)

		# store output
		self.codon_alignment = codon_alignment
	
	def write_fa(self) -> None:
		seqids = [self.terminals[key] for key in self.terminals]
		seqs = [''.join(next(iter(s)) for s in self.sequences[seqid]) for seqid in seqids]
		if len(seqids) != len(seqs):
			raise ValueError("IDs and sequences must have the same length.")
		
		sys.stdout.write("Alignment:\n")
		with open(self.out, "w", encoding="utf-8") as f:
			for seq_id, seq in zip(seqids, seqs):
				clean_seq = seq.strip().upper()
				f.write(f">{seq_id}\n{clean_seq}\n")
				sys.stdout.write(f">{seq_id}\n{clean_seq}\n")
		sys.stdout.write(f"Alignment saved into {self.out}\n//\n")


def main():
	"""
	Utilizes functions from other programs detect ORFs, translate sequences
	to amino acids, perform pairwise alignment, and back-translate alignments and
	generate codon position statistics and visualizations.
	"""
	args = parse_args()
	msaligner = MSAligner(
		fasta_file=args.input,
		aln_file=args.output,
		ref=args.reference,
		match=args.match,
		mismatch=args.mismatch,
		indel=args.indel,
		ttable=args.table,
		ksize=args.kmer
	)
	history: List[str] = []
	msaligner.detect_orfs()
	msaligner.find_kmer_translate()
	msaligner.count_kmers()
	msaligner.build_kmer_matrix()
	msaligner.plot_kmer_heatmap()

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

	# Write final alignment to disk
	msaligner.write_fa()
	
	
if __name__ == "__main__":
	main()