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



    This test performs the alignment of test_seq1.fa with test_seq2.fa using BLOSUM62 and a gap open penalty of -10 
    with ane extension penalty of -1.
    It then tests that the 3 alignment matrices were correctly filled by asserting each constructed alignment matrices are equal (I'm actually using np.allclose, which returns 
    True if the element-wise comparison of values between the two matrices are similar within a very tiny threshold) to the corresponding expected alignment matrix I calculated by hand
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    alignment = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",-10,-1)
    alignment.align(seq1,seq2)
    
    expected_align_matrix = [  #alignment matrix I calculated by hand
       [ 0., -np.inf, -np.inf, -np.inf, -np.inf],
       [-np.inf,   5., -12., -12., -14.],
       [-np.inf, -11.,   4.,  -1.,  -6.],
       [-np.inf, -13.,  -8.,   5.,   4.]
       ]
    expected_gapA_matrix = [ #gapA matrix I calculated by hand
       [-10., -11., -12., -13., -14.],
       [-np.inf, -22.,  -6.,  -7.,  -8.],
       [-np.inf, -23., -17.,  -7.,  -8.],
       [-np.inf, -24., -18., -18.,  -6.]]
    expected_gapB_matrix = [ #gapB matrix I calculated by hand
       [-10., -np.inf, -np.inf, -np.inf, -np.inf],
       [-11., -22., -23., -24., -25.],
       [-12.,  -6., -17., -18., -19.],
       [-13.,  -7.,  -7., -12., -17.]
       ]

    #assert each of the alignment matrices I calculated by hand are equal (within a very small threshold) of the returned alignment matrix by the module.
    assert np.allclose(alignment._align_matrix,expected_align_matrix)
    assert np.allclose(alignment._gapA_matrix,expected_gapA_matrix)
    assert np.allclose(alignment._gapB_matrix,expected_gapB_matrix)

    check_backtrace_equals_score(alignment) #check the alignment score returned by the Class object equals the alignment score recalculated from the final seqA and seqB alignments

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

    check_backtrace_equals_score(alignment) #check the alignment score returned by the Class object equals the alignment score recalculated from the final seqA and seqB alignments







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

    check_backtrace_equals_score(alignment) #check the alignment score returned by the Class object equals the alignment score recalculated from the final seqA and seqB alignments

def check_backtrace_equals_score(alignment):
    """
    This function takes the final seqA and seqB alignments and uses them to re-calculate the alignment score. It then checks that the re-calculated alignment score is the same as the 
    alignment score returned from the Class object. his ensures backtracking and alignment agree with each other; if they don't agree with each other, that would be a sign of improper implementation
    
    Parameters:
        NeedlemanWunsch.align object, which conatins the seqA and seqB alignmenst as well as the alignment scores

    """

    def add_gap_penalty(alignment_seq,i,score):
    """
    Checks if a gap is a gap opening + extension, or just a gap extension and applies the appropriate penalty
    """
        if i >= 1:
            
            if alignment_seq[i-1] == '-': #if i>=1, it is possible a gap is extending rather than being just opened, so check the previous character for a gap
                score += alignment.gap_extend #if there was a gap at the previous position of the alignment, add extension penalty only
            else:
                score += (alignment.gap_open + alignment.gap_extend)
        else:
                score += (alignment.gap_open + alignment.gap_extend) #if i == 0 the gap must have just opened, so add gap open +
        return score




    score = 0
    for i in range(len(alignment.seqB_align)): #iterate through length of one of the alignments. Both are the same length
        seqA_align_char = alignment.seqA_align[i]
        seqB_align_char = alignment.seqB_align[i]
        
        if seqA_align_char != '-' and seqB_align_char != '-':
            score += alignment.sub_dict[(seqA_align_char,seqB_align_char)] #if there is no gap at this position, add match value to the alignment score
        
        elif seqA_align_char == '-':
            score = add_gap_penalty(alignment.seqA_align,i,score) #if there is a gap on the seqA alignment, check if is an extension or opening + extension and apply penalty accordingly
            
        else: #if seqB_align_char == '_'
            score = add_gap_penalty(alignment.seqB_align,i,score)#if there is a gap on the seqB alignment, check if is an extension or opening + extension and apply penalty accordingly

    assert score == alignment.alignment_score #test the score computed here from the final alignment == the alignment score returned by the Class
