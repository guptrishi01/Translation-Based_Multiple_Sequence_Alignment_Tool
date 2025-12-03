#!/usr/bin/env python3

"""
This module provides functions for simulating the Needleman Wunsch Algorithm
for two sequences.

It includes functions for initializing the matrix, filling it out via Needleman Wunsch,
and tracing it back to find the optimal alignment. The goal is align each pairing of sequences
for back-translation and statistic + visualizations.

Functions:
	- init_mat(seq1, seq2) - Initializes matrix for two sequences
	- fill_matrix(matrix, seq1_array, seq2_array, match, mismatch, indel) - Fills out matrix for two sequences based on Needleman-Wunsch algorithm
	- trace_matrix(matrix, seq1_array, seq2_array) - Finds optimal alignment of two sequences from filled out matrix
"""

import numpy as np
import sys

def init_mat(seq1: np.ndarray, seq2: np.ndarray) -> np.ndarray:
	"""
	Initialize the Needleman-Wunsch scoring matrix

	Args:
		seq1 (np.ndarray) : NumPy array of characters representing sequence 1
		seq2 (np.ndarray) : NumPy array of characters representing sequence 2

	Returns:
		(np.ndarray) : Zero-initialized matrix for filling and tracing
	"""
	
	if not isinstance(seq1, np.ndarray) or not isinstance(seq2, np.ndarray):
		sys.stderr.write("Sequences for alignment must be numpy arrays")
		raise ValueError("Sequences for alignment must be numpy arrays")

	rows = len(seq2) + 1
	cols = len(seq1) + 1

	if rows <= 1 or cols <= 1:
		sys.stderr.write("Sequences must contain at least one element.")
		raise ValueError("Sequences must contain at least one element.")

	return np.zeros((rows, cols), dtype=int)
	

def fill_matrix(
	matrix: np.ndarray,
	seq1_array: np.ndarray,
	seq2_array: np.ndarray,
	match: int = 1,
	mismatch: int = -1,
	indel: int = -1
) -> np.ndarray:
	"""
	Fill the scoring matrix using the Needleman–Wunsch down-pass algorithm.

	Args:
		matrix (np.ndarray) : Zero-initialized matrix
		seq1_array (np.ndarray) : NumPy array of characters representing sequence 1.
		seq2_array (np.ndarray) : NumPy array of characters representing sequence 2.
		match (int) : Score for a character match (default  is +1).
		mismatch : (int) : Score for a character mismatch (default is -1).
		indel (int) : Penalty for an insertion/deletion (default is -1).

	Returns:
		(np.ndarray) : Scoring matrix filled via Needleman-Wunsch

	"""
	if not all(isinstance(x, int) for x in (match, mismatch, indel)):
		sys.stderr.write("match, mismatch, and indel must be integers.")
		raise TypeError("match, mismatch, and indel must be integers")
	
	rows, cols = matrix.shape
	if matrix.shape != (len(seq2_array) + 1, len(seq1_array) + 1):
		sys.stderr.write("Matrix shape does not match sequence lengths")
		raise ValueError("Matrix shape does not match sequence lengths")


	# Initialize the first row and first column with indel penalties
	for i in range(1, rows):
		matrix[i, 0] = matrix[i - 1, 0] + indel
	for j in range(1, cols):
		matrix[0, j] = matrix[0, j - 1] + indel

	# Fill the matrix
	for i in range(1, rows):
		for j in range(1, cols):
			# Assign scores based on if characters match or not
			if seq2_array[i - 1] == seq1_array[j - 1]:
				score_diagonal = matrix[i - 1, j - 1] + match
			else:
				score_diagonal = matrix[i - 1, j - 1] + mismatch

			# Indel scores
			score_up = matrix[i - 1, j] + indel
			score_left = matrix[i, j - 1] + indel

			# Take maximum
			max_value = max(score_diagonal, score_up, score_left)
			matrix[i, j] = max_value

	return matrix

def trace_matrix(
	matrix: np.ndarray,
	seq1_array: np.ndarray,
	seq2_array: np.ndarray
) -> np.ndarray:
	"""
	Trace back through the filled Needleman–Wunsch matrix to reconstruct
	the optimal alignment path.

	Args:
		matrix (np.ndarray) : Filled Needleman–Wunsch scoring matrix.
		seq1_array (np.ndarray) : NumPy array of characters representing sequence 1.
		seq2_array (np.ndarray) : NumPy array of characters representing sequence 2.

	Returns:
		(np.ndarray) : A NumPy array with alignments of sequence 1 and 2:

	"""
	if matrix.size == 0:
		sys.stderr.write("Cannot trace back an empty matrix")
		raise ValueError("Cannot trace back an empty matrix")
	i, j = matrix.shape[0] - 1, matrix.shape[1] - 1
	aligned_seq1 = []
	aligned_seq2 = []

	while i > 0 or j > 0:
		# Handle boundaries first
		if i == 0:
			aligned_seq1.append(seq1_array[j - 1])
			aligned_seq2.append(set(['-']))
			j -= 1
			continue
		elif j == 0:
			aligned_seq1.append(set(['-']))
			aligned_seq2.append(seq2_array[i - 1])
			i -= 1
			continue

		current = matrix[i, j]
		diag = matrix[i - 1, j - 1]
		up = matrix[i - 1, j]
		left = matrix[i, j - 1]

		# Determine which move leads to current
		if current == diag + (1 if seq1_array[j - 1] == seq2_array[i - 1] else -1):
			# Diagonal move (priority 1)
			aligned_seq1.append(seq1_array[j - 1])
			aligned_seq2.append(seq2_array[i - 1])
			i -= 1
			j -= 1
		elif current == left - 1:
			# Move left (gap in seq2, priority 2)
			aligned_seq1.append(seq1_array[j - 1])
			aligned_seq2.append(set(['-']))
			j -= 1
		else:
			# Move up (gap in seq1, priority 3)
			aligned_seq1.append(set(['-']))
			aligned_seq2.append(seq2_array[i - 1])
			i -= 1

	aligned_seq1.reverse()
	aligned_seq2.reverse()

	return np.array([aligned_seq1, aligned_seq2])

