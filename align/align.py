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


        Parameters:
            seqA: str
                First string of characters to be aligned against seqB
            seqB: str
                Second string of characters to be aligned against seqA

        Returns:
            Tuple with alignment score, seqA alignment, and seqB alignment

        ******************************************************************************************************************
        ********** For all comments in this package M[i,j] notation refers to row i and column j of matrix M *************
        ******************************************************************************************************************



        Brief overview of what this method is doing and how backtracing is implemented (more in depth explanation below):
        1) Initialize alignment and gap matrices
        2) Fill in each entry of the gap and alignment matrix starting at the top left, moving left to right and top to bottom in each matrix. 
        3) When each entry is filled, a decision is made where the value (plus gap penalty or match score additions )in an entry in a matrix will come from the maximum value of optimal choices that were made at an earlier matrix-filling step. 
        4) When this decision is made, the name of the matrix from which the current score partially derives from, as well as the location within that matrix, is saved in the backtracing matrix that corresponds to the matrix the entry was saved in
        5) for example, if self._align_matrix[3,3] came from match_score(seqA[2],seqB[2]) + self.align_matrix[2,2] -- meaning the [2,2] entry in self._align_matrix was greater than the values at the same position in either of the gap matrices --
        then at the [3,3] position of self._back (the backtracing matrix corresponding to self._align_matrix) the tuple ("Align Matrix",2,2) will be saved. This matrix says that the value in the [3,3] position of self._align_matrix was calculated in part
        by taking the value at the [2,2] position of the "Align Matrix" which codes for self._align_matrix. 
        6) Suppose the max alignment score was in the bottom right entry of the self._gapA_matrix. To find the alignment, one could index the bottom right entry of self._back_A. This would point to the position and matrix from which the previous optimal decision was made
        The alignment at that position would be saved in self.seqA_align and self.seqB_align. Then the next backtracing step would be found at the same position of the previous optimal decision's backtracing matrix. And this would repeat until position [0,0] is reached

        What the align method is doing:
            A. Initialize Alignment matrix, Gap Matrices, and backtacing matrices with the initialize() function
                i) Set [0,0] position of _align_matrix as a 0, leave the rest of the values in this matrix as -inf
                ii) Set each element in the 0th row of _gapA_matrix equal to the self.gap_open penalty + the self.gap_extend penalty (multiplied by the extension length at that element)
                    a) This row represents the score of best alignment between 0 characters of A and j characters of B that ends in a gap in A, where j is a column.
                    b) Therefore, a gap is being aligned to j characters of B, and the alignment score is equal to the self.gap_open + j*self.gap_extend penalty
                    c) Leave the 0th column of the _gapA_matrix equal to -inf as this column represents the score of best alignment between i characters of A and 0 characters of B that ends in a gap in A, where i is a row
                    d) That is nonsensical, so leave those values as -inf
                iii) Likewise, set each element in the 0th column of _gapB_matrix equal to the self.gap_open penalty + the self.gap_extend penalty (multiplied by the extension length at that element)
                iv) While initializing the 0th row of _gapA_matrix and 0th column of _gapB_matrix, also initialize the same positions in the _back_A and _back_B backtracing matrices
                    a) All scores in the 0th row of _gapA_matrix come from a gap extension penalty that increases from left to right in the 0th row of _gapA_matrix. Therefore the 0th row (not including [0,0]) of the self._back_A matrix
                    contains a tuple in each entry of the form ("gapA Matrix", 0, j-1), where j is a column. For example, self._back_A[0,5] == ('gapA Matrix', 0, 4). This means that the entry at self._back_A[0,5] comes from 
                    the "gapA matrix" -- which is a string encoding for self._gapA_matrix -- at position [0,4]. Likewise, self._back_A[0,4] == ('gapA Matrix', 0, 3) because  self._gapA_matrix[0,4] == self._gapA_matrix[0,3] + self.gap_extend
                    b) Likewise, the 0th column (except for 0,0) self._back_B matrix is initialized to contain tuples of the form ("gapB Matrix",i-1,0) where i is a column. this means that each entry in the 0th column of self.gapB_matrix comes from the previous row in the same column of that matrix (plus a gap extension penalty)
                    c) the [0,0] entry is left empty (None) as backtracing recursion will stop once this position is reached
            B. Via a nested for loop, iterate through every entry after the 1st column and row of the matrices. At each position, fill the alignment matrix, two gap matrices, and three backtracing matrices using the fill_entry() function
                i) First, the function saves the pair of characters this entry corresponds to in the self._align_matrix as a tuple `current_pair` that can be used to index self.sub_dict for a match score
                ii) Then, the function creates a dictionary `choices_dict`. Every time an entry in the alignment or gap matrices is filled, an optimal decision has to be made where the scores at a previous step
                from the alignment and two gap matrices are considered and the maximum is selected. np.argmax will be used to select the max value by returning the index of that value. This index will be passed into `choices_dict`
                to save the name of the matrix used for the optimal decision to facilitate backtracing at a later step
                iii) The entry at the current position of the loop in the self._align_matrix and self._back are filled
                    iiia) The maximum value between the alignment matrix and two gap matrices one position diagonally to the upper left is selected using np.argmax. The index of the max value is returned. 
                    iiib) if there is a tie, the first element is returned. This means ties will favor alignment, then gapA, and gapB last due to the order they were passed into np.argmax
                    iiic) the index of the max value is used to select the max value, add it to the match score at the current position (using self.sub_dict) and input it into the current position of the self.align_matrix
                    iiid)the index of the max value is also used to select the matrix from which that value is partially derived via choices_dict. The current position in self._back is then filled with a tuple of the form:
                    (Name of the matrix used to form the optimal decision at the current position in self._align_matrix, row # of optimal decision, col # of optimal decision ). When backtracing, this will be used to point
                    to the location of the matrix used to help calculate the score going into the current entry of self._align_matrix and therefore facilitating backtracing the alignment
                iv) Likewise the entry at the current position of the loop is filled in self._gapB_matrix and self._back_B and then self._gapA_matrix and self._backA
                    iva) This procedure is extremely similar as (iii) except the positions considered for the previous step in the alignment and gap matrices differ. 
                    Additionally a gap extension penalty (and possibly a gap open penalty) are added
            C. Once the alignment matrix, two gap matrices, and the backtracing matrices are filled, backtracing commences. For details on this procedure, see the _backtrace() docstring

        

        
        
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

        
        def initialize(): #initialize Alignment matrix, Gap matrices, and backtracing matrices
            self._align_matrix[0,0] = 0 #set 0,0 position of align_matrix to 0


            for j in range(len(seqB)+1): # all matrices have dimension (len(seqA)+1, len(seqB)+1); use len(seqB)+1 to iterate through columns and initialize 0th row of self_gapA_matrix
                #initialize the jth col in row 0 of gapA_matrix equal to gap open plus length of extension penalty
                self._gapA_matrix[0,j] = self.gap_open + (j * self.gap_extend)
                if j > 0: #add pointers moving to the left of gapA Matrix's first row. Leave [0,0] empty so you stop backtracing once you get to [0,0]
                    self._back_A[0,j] = ('gapA Matrix', 0, j - 1) #The jth entry in the 0th row of the self._gapA_matrix came from the j-1th entry in the 0th row of self._gap_A_matrix (plus a gap extension penalty). 
                    #self._back_A therefore includes this information at position [0,j] via a tuple (Matrix previous optimal decision comes from, row position of that optimal decision, col position of that optimal decision)
                
            for i in range(len(seqA)+1): # all matrices have dimension (len(seqA)+1, len(seqB)+1); use len(seqB)+1 to iterate through rows and initialize 0th col of self_gapB_matrix
                #initialize the ith row in col 0 of gapB_matrix equal to gap open plus length of extension penalty
                self._gapB_matrix[i,0] = self.gap_open + (i * self.gap_extend)
                if i > 0: #add pointers moving to the left of gapB Matrix's first row. Leave [0,0] empty so you stop backtracing once you get to [0,0]
                    self._back_B[i,0] = ('gapB Matrix', i-1,0) #The ith entry in the 0th col of the self._gapB_matrix came from the i-1th entry in the 0th col of self._gap_B_matrix (plus a gap extension penalty). 
                    #self._back_B therefore includes this information at position [i,0] via a tuple (Matrix previous optimal decision comes from, row position of that optimal decision, col position of that optimal decision)
                
            
        def fill_entry(current_row,current_column):
            current_pair = seqA[current_row-1],seqB[current_column-1] #Actual alignment pair is 1 position off because of the 0th row and 0th column of matrix do not correspond to a character
            
            #The entries in this choices_dict are in the same order as will be supplied to np.argmax to select an optimal decision. 
            #This dictionary will therefore be used to convert the index returned from np.argmax to a string with the name of the corresponding matrix holding the optimal decision in the previous step
            #This string holding the name of the matrix will be used for backtracing
            self.choices_dict  = {  
                0: 'Align Matrix',
                1: 'gapA Matrix',
                2: 'gapB Matrix'
            }

            ############################## Fill entry for self._align_matrix and self._back #####################################################
            align_matrix_decision_choices = [
                self._align_matrix[current_row - 1,current_column - 1], # consider score from align_matrix at i-1, j-1 position
                self._gapA_matrix[current_row - 1,current_column - 1], #consider score from gapA_matrix at i-1,j-1 position
                self._gapB_matrix[current_row - 1,current_column - 1] #consider score from gapB_matrix at i-1,j-1 position
                ]

            optimal_decision_align_matrix_idx =  np.argmax(align_matrix_decision_choices) #get index corresponding to choice with max value. If two options have same value, pick the first one, which means the align matrix will be favored in a tie, 
            #followed by the gapA matrix, with the gapB matrix coming last in a tie.
            
            #value at this current entry in align_matrix will be equal to the match score of the current_pair + max of the 3 options
            self._align_matrix[current_row,current_column] = self.sub_dict[current_pair] + align_matrix_decision_choices[optimal_decision_align_matrix_idx]

            # At the position of the current entry, add to the self._back pointer matrix a tuple of form:
            # (Name of the matrix used to form the optimal decision at the current position in self._align_matrix, row # of optimal decision, col # of optimal decision )
            # Name of matrix used to form the optimal decision comes from choices_dict. The index chosen by np.argmax corresponds to this matrix, and the dictionary converts that index into a string            
            self._back[current_row,current_column] = (self.choices_dict[optimal_decision_align_matrix_idx],current_row - 1,current_column - 1)
            


            ############################## Fill entry for self._gapB_matrix and self._back_B #####################################################
            gapB_matrix_decision_choices = [
                    self.gap_open + self.gap_extend + self._align_matrix[current_row - 1,current_column], #consider score from previous row, same col of align matrix, plus gap open and extend penalty
                    self.gap_open + self.gap_extend + self._gapA_matrix[current_row - 1,current_column], #consider score from previous row, same col of gapA_matrix, plus gap open and extend penalty
                    self.gap_extend + self._gapB_matrix[current_row - 1,current_column] #consider score from previous row, same col of align matrix, plus gap extend. The gap is already open and continuing in this case, so no gap open penalty
                    ]
            optimal_decision_gapB_matrix_idx = np.argmax(gapB_matrix_decision_choices) #get index of optimal score. Ties will go to element with highest index. 
            
            
            self._gapB_matrix[current_row,current_column] = gapB_matrix_decision_choices[optimal_decision_gapB_matrix_idx] #add optimal score to gapB_matrix at current position
            self._back_B[current_row,current_column] = (self.choices_dict[optimal_decision_gapB_matrix_idx],current_row - 1,current_column) #add tuple with (name of matrix used to form optimal decision, row# holding entry for optimal decision, col# holding entry for optimal decision) to self._back_B
           


            ############################## Fill entry for self._gapA_matrix and self._back_A #####################################################
            gapA_matrix_decision_choices = [
                        self.gap_open + self.gap_extend + self._align_matrix[current_row ,current_column -1], #consider score from previous col, same row of align matrix, plus gap open and extend penalty
                        self.gap_extend + self._gapA_matrix[current_row, current_column - 1], #consider score from previous col, same row of align matrix, plus extend penalty. Gap is already open here so no gap open penalty
                        self.gap_open + self.gap_extend + self._gapB_matrix[current_row, current_column - 1] #consider score from previous col, same row of gapB matrix, plus gap open and extend penalty
                        ]
            optimal_decision_gapA_matrix_idx  = np.argmax(gapA_matrix_decision_choices)
            
            self._gapA_matrix[current_row,current_column] = gapA_matrix_decision_choices[optimal_decision_gapA_matrix_idx] #add optimal score to gapA_matrix at current position
            self._back_A[current_row,current_column] = (self.choices_dict[optimal_decision_gapA_matrix_idx], current_row, current_column - 1) #add tuple with (name of matrix used to form optimal decision, row# holding entry for optimal decision, col# holding entry for optimal decision) to self._back_A
           





                
        initialize() # call initialize function to initialize gap and alignment matrices
        
        #fill each entry in all 6 matrices
        for row in range(len(seqA) + 1):
            for column in range(len(seqB) + 1):
                if row !=0 and column != 0: #only fill entry outside of 0th row and 0 column; these are already initialized
                    fill_entry(row,column) #fills alignment matrix, gap matrix, and their respective backtracing pointer matrices
        
        
        
        
                                           
                                           
       

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        # TODO Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.



        Returns:
            Tuple with alignment score, seqA alignment, and seqB alignment: 




        What is happening in the _backtrace method:
            A. Take the alignment score in the bottom right of the 3 matrices holding the alignment scores in a list  `final_alignment_score_choices`
                i) the bottom right of these matrices holds the cumulative alignment score for each matrix. So the max of these three options gives the max total alignment score, and backtracing 
                through that matrix starting at the bottom right will yield the optimal alignment
            B) use np.argmax to find the index of that list with the max total alignment score. If there is a tie, the self._align_matrix will be favored, followed by self._gapA_matrix, followed by self._gapB_matrix, 
                because np.argmax returns the first index in the event of a tie
            C) Take the index with the max value and use it to index `final_alignment_score_choices` and save the max alignment score in self.alignment_score
            D) Begin recrusively_backtrace() to backtrace the alignment
                i) Call the function by inputing the name of the matrix that had the max alignment score. Use self.choices_dict to index the name of this matrix; the keys of this dictionary line up with the matrices in 
                `final_alignment_score_choices`. Also include the bottom right row,column position as the second and third arguments
                ii) The base case for this recursive function is when you reach the [0,0] entry during backtracing. Once you arrive here, recursion will stop and backtracing will be complete
                iii) An if/elif statement splits this function based on the name of the matrix through which backtracing is commencing. 
                    iiia) For example, if the function was called using the "Align Matrix" at position [4,4], the characters in seqA and seqB corresponding to the current position of the matrix  are added to the front of the alignment, as backtracing is currently at the alignment matrix
                        iiiaa) Because there is a 0 row and 0 column in the matrix, position [4,4] corresponds to the character at index 3 of seqA and index 3 of seqB. The characters at this position are added to the front of self._seqA and self._seqB, respectively
                    iiib) Additionally, the current position (position [4,4]) of the corresponding backtracing matrix, self._back, is accessed, returning a tuple that points to the position of the optimal decision in the previous step
                        iiiba) For example, this could return ("gapA Matrix",3,3), which would mean the alignment score at self._align_matrix[4,4] came from adding the match score at [4,4] to the cumulative score in self._gapA_matrix[3,3]
                    iiic) Call recursively_backtrace again using the entries of this tuple as arguments to initiate recursion. If, as in the example, ("gapA Matrix",3,3) was used then a gap in seqA is added to the front of the alignment 
                    and the corresponding character in seqB is added to the seqB alignment. Additionally, the self._back_A backtracing matrix is accessed at position [3,3] to get the matrix and position from which the previous optimal score came
                    iiid) This repeats until the [0,0] position is reached, terminating recursion. The alingment for seqA and seqB will have been stored in seqA_align and seqB_align.

        """        
            
        # Put alignment scores in the bottom right of the three alignment matrices in a list. One of these is the max score for the total alignment. 
        final_alignment_score_choices = [
                    self._align_matrix[len(self._seqA),len(self._seqB)],
                    self._gapA_matrix[len(self._seqA),len(self._seqB)],
                    self._gapB_matrix[len(self._seqA),len(self._seqB)],
                ]

        
        final_alignment_optimal_score_idx = np.argmax(final_alignment_score_choices) #get index corresponding to max alignment score for total alignment
        # If there is a tie, np.argmax will favor the alignment from align_matrix, then gapA, then gapB


        self.alignment_score = final_alignment_score_choices[final_alignment_optimal_score_idx] #Save max alignment score in attribute
        
        def recursively_backtrace(score_matrix_name,score_matrix_row,score_matrix_column):
            if score_matrix_row != 0 or score_matrix_column != 0: #Base case. Stop recursing once you reach top left of backtrace
                if score_matrix_name == "Align Matrix":
                    self.seqA_align = self._seqA[score_matrix_row - 1] + self.seqA_align # add matching element of seqA to front of seqA alignment that is being recursively built
                    self.seqB_align = self._seqB[score_matrix_column - 1] + self.seqB_align # treat seqB accordingly, since the optimal score concluded seqA and seqB matched/mismatched here
                    next_args = self._back[score_matrix_row,score_matrix_column]
                    return recursively_backtrace(next_args[0],next_args[1],next_args[2]) #Inputting the current row and column into self._back returns a tuple with
                #the score matrix name, row, and column, of the optimal decision at the previous step, facilitating backtrace
            
                elif score_matrix_name == "gapA Matrix":
                    self.seqA_align = '-' + self.seqA_align #add gap in A alignment
                    self.seqB_align = self._seqB[score_matrix_column - 1] + self.seqB_align #do not add gap in B alignment
                    next_args = self._back_A[score_matrix_row,score_matrix_column]

                    return recursively_backtrace(next_args[0],next_args[1],next_args[2]) #Inputting the current row and column into self._back_A returns a tuple with
                #the score matrix name, row, and column, of the optimal decision at the previous step, facilitating backtrace

                elif score_matrix_name == "gapB Matrix":
                    self.seqA_align = self._seqA[score_matrix_row - 1] + self.seqA_align #do not add gap in A alignment
                    self.seqB_align = '-' + self.seqB_align #add gap in B alignment
                    next_args = self._back_B[score_matrix_row,score_matrix_column]
                    return recursively_backtrace(next_args[0],next_args[1],next_args[2]) #Inputting the current row and column into self._back_B returns a tuple with
                #the score matrix name, row, and column, of the optimal decision at the previous step, facilitating backtrace
            
        recursively_backtrace(self.choices_dict[final_alignment_optimal_score_idx],len(self._seqA),len(self._seqB)) #begin backtracing at bottom right of score matrix that contains the max alignment score
        #The name of this matrix inputed into the first argument comes from accessing self.choices_dict using the `final_alignment_optimal_score_idx`. The order of the names in self.choices_dict is in the same order 
        #as they appeared in `final_alignment_score_choices`, so the index corresponding to matrix with the max alignment score was accessed in this way.
        
        return (self.alignment_score, self.seqA_align,self.seqB_align) #return tuple with the alignment score, seqA alignment, and seqB alignment
        

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
