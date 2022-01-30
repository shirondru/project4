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
    alignment = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",-10,-1)
    alignment.align(seq3,seq4)

    assert alignment.alignment_score == 17, "Your Alignment Score is Different than Expected" #assert alignment score is 17, the value the TAs said to expect
    assert alignment.seqA_align == "MAVHQLIRRP", "Seq3 Alignment is wrong" #assert Seq3 alignment is what the TAs said to expect
    assert alignment.seqB_align == "M---QLIRHP", "Seq4 Alignment is wrong" #assert Seq4 alignment is what the TAs said to expect

    






def test_example_from_class():
    """
    Asserting that I get same alignment score as the example shown in lecture
    Using a matrix `test_sub_mat` that I created that scores -1 for mismatches and 1 for correct matches, just like the lecture example
    Gap open = -3, gap extend = -1, just like the lecture example

    Importantly, this also tests my method favors matches/mismatches over gaps when there is a tie in alignment scores, as per the heuristic I introduced
    #this example could have also given "AA--T" for the seqY alignment due to ties in the alignment scores. But it should return "--AAT"
    """

    seqX = "ACACT" 
    seqY = "AAT"
    alignment = NeedlemanWunsch("./substitution_matrices/test_sub_mat.mat",-3,-1)
    alignment.align(seqX,seqY)
    assert alignment.alignment_score == -4, "Your alignment score is different than expected" #Class example had alignment score -4, assert you were able to reproduce that
    assert alignment.seqA_align == "ACACT", "seqX Alignment is wrong" #Class example had this alignment as ACACT, assert you were able to reproduce that
    assert alignment.seqB_align == "--AAT", "seqY Alignment is wrong" #Class example had this alignment as "--AAT", assert you were able to reproduce that
    #Additionally, this tests that ties in alignment are handled properly per my heuristic

def check_backtrace_equals_score():
    """
    Take the alignment output and use it to calculate a score. Check that score equals the alignment score output to ensure backtracing and alignment were implemented properly. If they don't corroborate each other, there is a problem
    """

    pass
    
