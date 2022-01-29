# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        # TODO: Fill in the Needleman-Wunsch Algorithm below
        to perform global sequence alignment of seqA and seqB
        and return a tuple with the following format
        (alignment score, seqA alignment, seqB alignment)
        Also, write up a docstring for this function using the
        _read_sub_matrix as an example.
        Don't forget to comment your code!
        
        
        The entry in align_matrix[i,j] is the maximum score for the alignment between seqA[1:i] and seqB[1:j]
        """
        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        self._back = np.empty((len(seqA) + 1, len(seqB) + 1),tuple)
        self._back_A = np.empty((len(seqA) + 1, len(seqB) + 1),tuple)
        self._back_B = np.empty((len(seqA) + 1, len(seqB) + 1),tuple)

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # TODO Implement the global sequence alignment here
        def initialize():
            self._align_matrix[0,0] = 0 #set 0,0 position of align_matrix to 0


            for j in range(len(seqB)+1):
                #initialize the jth col in row 0 of gapA_matrix equal to gap open plus length of extension penalty
                self._gapA_matrix[0,j] = self.gap_open + (j * self.gap_extend)
                if j > 0: #add pointers moving to the left of gapA Matrix's first row. Leave [0,0] empty so you stop backtracing once you get to [0,0]
                    self._back_A[0,j] = ('gapA Matrix', 0, j - 1)
                
            for i in range(len(seqA)+1):
                #initialize the ith row in col 0 of gapB_matrix equal to gap open plus length of extension penalty
                self._gapB_matrix[i,0] = self.gap_open + (i * self.gap_extend)
                if i > 0: #add pointers moving to the left of gapB Matrix's first row. Leave [0,0] empty so you stop backtracing once you get to [0,0]
                    self._back_B[i,0] = ('gapB Matrix', 0, i - 1)
                
            
        def fill_entry(current_row,current_column):
            current_pair = seqA[current_row-1],seqB[current_column-1] #Actual alignment pair is 1 position off because of the 0th row and 0th column of matrix do not correspond to a character
            
            
            #max built-in method will return max score at earlier step along with the matrix that score came from in a tuple of the form (max score, name of the matrix the score comes from)
            #The 0th index (i.e, the value) of the tuple with the max value will be indexed to fill into the current entry in the self._align_matrix
            #The 1th index (i.e, the name of the matrix the score came from) will be used to build the self._back pointer matrix
            optimal_decision_align_matrix =  max(
                (self._align_matrix[current_row - 1,current_column - 1],'Align Matrix'), # consider score from align_matrix at i-1, j-1 position
                (self._gapA_matrix[current_row - 1,current_column - 1],'gapA Matrix'), #consider score from gapA_matrix at i-1,j-1 position
                (self._gapB_matrix[current_row - 1,current_column - 1],'gapB Matrix') #consider score from gapB_matrix at i-1,j-1 position
                )
            
            #value at this current entry in align_matrix will be equal to the match score of the current_pair + max of the 3 options
            self._align_matrix[current_row,current_column] = self.sub_dict[current_pair] + optimal_decision_align_matrix[0]
            #At the position of the current entry, add to the self._back pointer matrix a tuple of form:
            # (Name of the matrix used to form the optimal decision at the current position in self._align_matrix, row # of optimal decision, col # of optimal decision )
            self._back[current_row,current_column] = (optimal_decision_align_matrix[1],current_row - 1,current_column - 1)
            
            optimal_decision_gapB_matrix = max(
                (self.gap_open + self.gap_extend + self._align_matrix[current_row - 1,current_column],'Align Matrix'),
                (self.gap_extend + self._gapB_matrix[current_row - 1,current_column],'gapB Matrix'),
                (self.gap_open + self.gap_extend + self._gapA_matrix[current_row - 1,current_column],'gapA Matrix')
            )
            
            
            self._gapB_matrix[current_row,current_column] = optimal_decision_gapB_matrix[0]
            self._back_B[current_row,current_column] = (optimal_decision_gapB_matrix[1],current_row - 1,current_column)
            
            
            optimal_decision_gapA_matrix = max(
                (self.gap_open + self.gap_extend + self._align_matrix[current_row ,current_column -1],'Align Matrix'),
                (self.gap_open + self.gap_extend + self._gapB_matrix[current_row, current_column - 1],'gapB Matrix'),
                (self.gap_extend + self._gapA_matrix[current_row, current_column - 1],'gapA Matrix'),
            )
            
            self._gapA_matrix[current_row,current_column] = optimal_decision_gapA_matrix[0]
            self._back_A[current_row,current_column] = (optimal_decision_gapA_matrix[1], current_row, current_column - 1)
                
        initialize() # call initialize function to initialize gap and alignment matrices
        
        #fill each entry in all 6 matrices
        for row in range(len(seqA) + 1):
            for column in range(len(seqB) + 1):
                if row !=0 and column != 0: #only fill entry outside of 0th row and 0 column
                    fill_entry(row,column) #fills alignment matrix, gap matrix, and their respective pointer matrices
        
        
        
        
                                           
                                           
       

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        # TODO Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.
        """
        # Implement this method based upon the heuristic chosen in the align method above.
        
            
        # Find optimal final score
        optimal_decision = max(
        (self._align_matrix[len(seqA),len(seqB)],'Align Matrix'),
        (self._gapA_matrix[len(seqA),len(seqB)],'gapA Matrix'),
        (self._gapB_matrix[len(seqA),len(seqB)],'gapB Matrix')  
        )
        
        self.alignment_score = optimal_decision[0] #max score from bottom right of all 3 matrices is the alignment score
        
        def recursively_backtrace(score_matrix_name,score_matrix_row,score_matrix_column):
            if score_matrix_row != 0 and score_matrix_column != 0: #Base case. Stop recursing once you reach top left of backtrace
                if score_matrix_name == "Align Matrix":
                    self.seqA_align = seqA[score_matrix_row - 1] + self.seqA_align # add matching element of seqA to front of seqA alignment that is being recursively built
                    self.seqB_align = seqB[score_matrix_column - 1] + self.seqB_align # treat seqB accordingly, since the optimal score concluded seqA and seqB matched/mismatched here
                    next_args = self._back[score_matrix_row,score_matrix_column]
                    return recursively_backtrace(next_args[0],next_args[1],next_args[2]) #Inputting the current row and column into self._back returns a tuple with
                #the score matrix name, row, and column, of the optimal decision at the previous step, facilitating backtrace
            
                elif score_matrix_name == "gapA Matrix":
                    self.seqA_align = '-' + self.seqA_align #add gap in A alignment
                    self.seqB_align = seqB[score_matrix_column - 1] + self.seqB_align #do not add gap in B alignment
                    next_args = self._back_A[score_matrix_row,score_matrix_column]

                    return recursively_backtrace(next_args[0],next_args[1],next_args[2]) #Inputting the current row and column into self._back_A returns a tuple with
                #the score matrix name, row, and column, of the optimal decision at the previous step, facilitating backtrace

                elif score_matrix_name == "gapB Matrix":
                    self.seqA_align = seqA[score_matrix_row - 1] + self.seqA_align #do not add gap in A alignment
                    self.seqB_align = '-' + self.seqB_align #add gap in B alignment
                    next_args = self._back_B[score_matrix_row,score_matrix_column]
                    return recursively_backtrace(next_args[0],next_args[1],next_args[2]) #Inputting the current row and column into self._back_B returns a tuple with
                #the score matrix name, row, and column, of the optimal decision at the previous step, facilitating backtrace
            
        recursively_backtrace(optimal_decision[1],len(seqA),len(seqB)) #begin backtracing at bottom right of score matrix that contains the max alignment score
        
        
        #1. align_matrix has max score in bottom right. call function with "Align Matrix",4,6
        #2. First if statement is true. set T and T for alignment of both seqA and seqB
        #3 lookup this position in self._back. This returns the next matrix name and row and column
def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
