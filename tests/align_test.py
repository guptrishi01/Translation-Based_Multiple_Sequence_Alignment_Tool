#!/usr/bin/env python3

import pytest
import pandas as pd
import numpy as np
from msaligner.needleman_wunsch import init_mat, fill_matrix, trace_matrix

def test_init_mat_basic():
    """
    Testing initialization of matrix
    """
    s1 = np.array(list("ABC"))
    s2 = np.array(list("XYZ"))
    m = init_mat(s1, s2)
    assert m.shape == (4, 4)
    assert np.all(m == 0)

def test_fill_matrix_basic():
    """
    Testing the correctness of filling a matrix
    """
    s1 = np.array(list("GATT"))
    s2 = np.array(list("GCT"))
    m = init_mat(s1, s2)
    m2 = fill_matrix(m, s1, s2)

    assert m2.shape == (len(s2)+1, len(s1)+1)
    # Check that scoring was applied
    assert np.any(m2[1:, 1:] != 0)


def test_trace_simple_alignment():
    """
    Testing the tracing of matrix
    """
    s1 = np.array(list("GATTACA"))
    s2 = np.array(list("GCATGCU"))

    m = fill_matrix(init_mat(s1, s2), s1, s2)
    aln = trace_matrix(m, s1, s2)

    assert aln.shape == (2, len(aln[0]))
    assert len(aln[0]) == len(aln[1])

def test_trace_matrix_correctness():
    """
    Testing correctness of tracing a matrix
    """
    s1 = np.array(list("AA"))
    s2 = np.array(list("A"))

    m = fill_matrix(init_mat(s1, s2), s1, s2)
    aln = trace_matrix(m, s1, s2)

    s1 = "".join(next(iter(x)) for x in aln[0])
    s2 = "".join(next(iter(x)) for x in aln[1])

    # Should align A with A and introduce one gap
    assert s1 in ["AA"]
    assert s2.count("-") == 1
    assert s2.replace("-", "") == "A"

