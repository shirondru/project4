# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    pass

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    alignment = NeedlemanWunsch("./substitution_matrices/BLOSUM62",-10,-1)
    alignment.align(seq3,seq4)

    assert alignment.alignment_score == 17, "Your Alignment Score is Different than Expected" #assert alignment score is 17, the value the TAs said to expect
    assert alignment.seqA_align == "MAVHQLIRRP", "Seq3 Alignment is wrong" #assert Seq3 alignment is what the TAs said to expect
    assert alignment.seqB_align == "M---QLIRHP", "Seq4 Alignment is wrong" #assert Seq4 alignment is what the TAs said to expect

    assert






def test_example_from_class():
    """
    Asserting that I get same alignment score as the example shown in lecture
    """




